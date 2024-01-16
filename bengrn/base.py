"""
bengrn base module.
"""

import pandas as pd
import urllib.request
import os.path

from grnndata import GRNAnnData, from_adata_and_longform, from_scope_loomfile, utils
from anndata import AnnData, concat
from typing import Optional
from .tools import GENIE3
import numpy as np

from arboreto.algo import grnboost2

from pyscenic.utils import modules_from_adjacencies
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell

from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase

import seaborn as sns
import decoupler as dc
import scipy.stats
from sklearn.metrics import precision_recall_curve, auc, PrecisionRecallDisplay


# example constant variable
NAME = "bengrn"

FILEDIR = os.path.dirname(os.path.realpath(__file__))


class BenGRN:
    def __init__(self, grn: GRNAnnData, full_dataset: Optional[AnnData] = None):
        """
        Initializes the BenGRN class.

        Parameters:
            grn (GRNAnnData): The Gene Regulatory Network data.
            full_dataset (Optional[AnnData]): The full dataset, defaults to None.
        """
        self.grn = grn
        self.full_dataset = full_dataset

    def do_tests(self, to: str = "collectri", organism="human"):
        """
        This method performs tests on the Gene Regulatory Network (GRN) and the full dataset.
        It compares the GRN network structure with the ground truth (GT) database and computes various metrics.

        Parameters:
            to (str, optional): The name of the GT database to use. Defaults to "collectri".
            organism (str, optional): The organism to consider for the GT database. Defaults to "human".

        Returns:
            dict: A dictionary containing the computed metrics.
        """
        self.gt = get_GT_db(name=to, organism=organism)
        # similarity in network structure
        # Get the GRN network from both grn objects
        grn_self = self.grn.grn
        grn_gt = self.gt.grn

        intersection = grn_self.index.intersection(grn_gt.index)
        print("intersection of {} genes".format(intersection.size))
        res = {"intersection": intersection.size / grn_self.index.size}
        intersection = list(intersection)
        grn_self = grn_self.loc[intersection, intersection]
        grn_gt = grn_gt.loc[intersection, intersection]
        similar_edges = (
            (grn_self != 0) & (grn_gt != 0) & (np.sign(grn_self) == np.sign(grn_gt))
        )
        similar_edges_ct = np.sum(np.sum(similar_edges))

        # Compute the total number of edges
        total_edges = np.sum(np.sum(grn_self != 0))

        # Compute the total number of edges in the other GRN
        total_edges_other = np.sum(np.sum((grn_gt != 0)))
        precision = similar_edges_ct / total_edges
        recall = similar_edges_ct / total_edges_other
        accuracy = similar_edges_ct / (
            total_edges + total_edges_other - similar_edges_ct
        )

        # Compute the Spearman's rank correlation between the two overlapping sets of edges
        spearman_corr = scipy.stats.spearmanr(
            grn_self.values[similar_edges].flatten(),
            grn_gt.values[similar_edges].flatten(),
        )
        res.update(
            {
                "precision": precision,
                "recall": recall,
                "accuracy": accuracy,
                "spearman_corr": spearman_corr,
            }
        )

        gt = grn_gt.values.flatten() != 0
        preval = gt.sum() / len(gt)
        pr, re, _ = precision_recall_curve(gt, grn_self.values.flatten())
        res.update({"auc": auc(re, pr)})
        disp = PrecisionRecallDisplay(pr, re, prevalence_pos_label=preval)
        disp.plot(plot_chance_level=True)
        disp.ax_.set_ylim([0, min(pr[:-1].max() * 3, 1)])

        return res

    def get_self_metrics(self):
        print("I'm getting my own metrics")
        pass

        # recap N protein complexes

        # gene pairs reg correlation with high string score vs none

        # recap well-known relationships

        #

    def compare_to(self, other: GRNAnnData):
        print("I'm comparing myself to another dataset")
        pass


def get_scenicplus(
    filepath=FILEDIR + "/../data/10xPBMC_homo_scenicplus_genebased_scope.loom",
):
    """
    This function retrieves a loomx scenicplus data from a given file path and loads it as a GrnnData

    Parameters:
        filepath : str, optional
            The path to the scenicplus data file.
            Default is FILEDIR + "/../data/10xPBMC_homo_scenicplus_genebased_scope.loom".

    Raises:
        FileNotFoundError: If the file at the given path does not exist.

    Returns:
        GrnAnnData: The scenicplus data from the given file as a grnndata object
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(
            f"The file {filepath} does not exist. You likely need to download \
                this or another loomxfile from the scope website"
        )

    return from_scope_loomfile(filepath)


def get_sroy_gt(join="outer"):
    """
    This function retrieves the ground truth data from Stone and Sroy's study.

    Parameters:
        join : str, optional
            The type of join to be performed when concatenating the data.
            Default is "outer".

    Returns:
        GrnAnnData: The ground truth data as a grnndata object
    """
    adata_liu = AnnData(
        (
            2
            ** pd.read_csv(
                FILEDIR
                + "/../data/GroundTruth/stone_and_sroy/scRNA/liu_rna_filtered_log2.tsv",
                sep="\t",
            )
        )
        - 1
    ).T
    adata_chen = AnnData(
        (
            2
            ** pd.read_csv(
                FILEDIR
                + "../data/GroundTruth/stone_and_sroy/scRNA/chen_rna_filtered_log2.tsv",
                sep="\t",
            )
        )
        - 1
    ).T
    df = pd.read_csv(
        FILEDIR + "/../data/GroundTruth/stone_and_sroy/hESC_ground_truth.tsv",
        sep="\t",
        header=None,
    )
    adata_liu.obs["dataset"] = "liu"
    adata_chen.obs["dataset"] = "chen"
    adata = concat([adata_liu, adata_chen], join=join)
    return from_adata_and_longform(adata, df)


def compute_scenic(adata, data_dir=FILEDIR + "/../data"):
    """
    This function computes the SCENIC algorithm on the given data.

    Parameters:
        adata (AnnData): The annotated data matrix of shape n_obs x n_vars. Rows correspond to cells and columns to genes.
        data_dir (str, optional): The directory where the data files will be stored. Defaults to FILEDIR + "/../data".

    Returns:
        GRNAnnData: The Gene Regulatory Network data.
    """
    os.makedirs(data_dir, exist_ok=True)

    url1 = "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
    url2 = "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"

    file1 = os.path.join(data_dir, os.path.basename(url1))
    file2 = os.path.join(data_dir, os.path.basename(url2))

    if not os.path.exists(file1):
        urllib.request.urlretrieve(url1, file1)
    if not os.path.exists(file2):
        urllib.request.urlretrieve(url2, file2)
    dbs = RankingDatabase(fname=file2, name="genes_10kb_human_rank_v10")

    adjacencies = grnboost2(
        expression_data=adata.to_df(), tf_names=utils.TF, verbose=True
    )
    modules = list(modules_from_adjacencies(adjacencies, adata.to_df()))

    df = prune2df([dbs], modules, file1)
    regulons = df2regulons(df)
    # auc_mtx = aucell(adata.to_df(), regulons, num_workers=8)
    # sns.clustermap(auc_mtx, figsize=(12,12))

    # compute the grn matrix
    var_names = adata.var_names.tolist()
    da = np.zeros((len(var_names), len(var_names)), dtype=np.float)
    for reg in regulons:
        i = var_names.index(reg.transcription_factor)
        for k, v in reg.gene2weight.items():
            da[i, var_names.index(k)] = v
    # create the grndata
    grn = GRNAnnData(adata, grn=da)

    da = np.zeros((len(var_names), len(var_names)), dtype=np.float)
    for reg in modules:
        i = var_names.index(reg.transcription_factor)
        for k, v in reg.gene2weight.items():
            da[i, var_names.index(k)] = v
    grn.varp["modules"] = da
    return grn


def compute_genie3(adata, nthreads=30, ntrees=100, **kwargs):
    """
    This function computes the GENIE3 algorithm on the given data.

    Parameters:
        adata (AnnData): The annotated data matrix of shape n_obs x n_vars. Rows correspond to cells and columns to genes.
        nthreads (int, optional): The number of threads to use for computation. Defaults to 30.
        ntrees (int, optional): The number of trees to use for the Random Forests. Defaults to 100.
        **kwargs: Additional arguments to pass to the GENIE3 function.

    Returns:
        GRNAnnData: The Gene Regulatory Network data computed using the GENIE3 algorithm.
    """
    mat = np.asarray(adata.X.todense())
    names = adata.var_names[mat.sum(0) > 0].tolist()
    var = adata.var[mat.sum(0) > 0]
    mat = mat[:, mat.sum(0) > 0]
    VIM = GENIE3(mat, gene_names=names, nthreads=nthreads, ntrees=ntrees, **kwargs)
    return GRNAnnData(grn=VIM, X=mat, var=var, obs=adata.obs)


def get_GT_db(name="collectri", cell_type=None, organism="human", split_complexes=True):
    """
    use_prior_network loads a prior GRN from a list of available networks.

    Args:
        name (str, optional): name of the network to load. Defaults to "collectri".
        organism (str, optional): organism to load the network for. Defaults to "human".
        split_complexes (bool, optional): whether to split complexes into individual genes. Defaults to True.

    Raises:
        ValueError: if the provided name is not amongst the available names.
    """
    # TODO: use omnipath instead
    if name == "tflink":
        TFLINK = "https://cdn.netbiol.org/tflink/download_files/TFLink_Homo_sapiens_interactions_All_simpleFormat_v1.0.tsv.gz"
        net = pd_load_cached(TFLINK)
        net = net.rename(columns={"Name.TF": "source", "Name.Target": "target"})
    elif name == "htftarget":
        HTFTARGET = "http://bioinfo.life.hust.edu.cn/static/hTFtarget/file_download/tf-target-infomation.txt"
        net = pd_load_cached(HTFTARGET)
        net = net.rename(columns={"TF": "source"})
    elif name == "collectri":
        net = dc.get_collectri(organism=organism, split_complexes=split_complexes).drop(
            columns=["PMID"]
        )
    elif name == "dorothea":
        net = dc.get_dorothea(organism=organism)
    elif name == "regulonDB":
        pass
    else:
        raise ValueError(f"provided name: '{name}' is not amongst the available names.")
    varnames = list(set(net.iloc[:, :2].values.flatten()))
    adata = AnnData(var=varnames)
    adata.var_names = varnames
    grn = from_adata_and_longform(adata, net, has_weight=True)
    return grn


def pd_load_cached(url, loc="/tmp/", cache=True, **kwargs):
    # Check if the file exists, if not, download it
    loc += url.split("/")[-1]
    if not os.path.isfile(loc) or not cache:
        urllib.request.urlretrieve(url, loc)
    # Load the data from the file
    return pd.read_csv(loc, **kwargs)
