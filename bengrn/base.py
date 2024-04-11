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

import decoupler as dc
import scipy.stats
from sklearn.metrics import precision_recall_curve, auc, PrecisionRecallDisplay
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split

from omnipath.interactions import AllInteractions
from omnipath.requests import Annotations

import gseapy as gp

import tqdm

import matplotlib.pyplot as plt

# example constant variable
NAME = "bengrn"

FILEDIR = os.path.dirname(os.path.realpath(__file__))


# Define a function to calculate precision and recall
def precision_recall(true_con, grn_con):
    tp = len(true_con & grn_con)
    fp = len(grn_con - true_con)
    fn = len(true_con - grn_con)

    precision = tp / (tp + fp)
    recall = tp / (tp + fn)

    return precision, recall


class BenGRN:
    def __init__(
        self, grn: GRNAnnData, full_dataset: Optional[AnnData] = None, doplot=True
    ):
        """
        Initializes the BenGRN class.

        Parameters:
            grn (GRNAnnData): The Gene Regulatory Network data.
            full_dataset (Optional[AnnData]): The full dataset, defaults to None.
        """
        self.grn = grn
        self.full_dataset = full_dataset
        self.doplot = doplot

    def do_tests(
        self,
    ):
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

    def compare_to(
        self,
        other: Optional[GRNAnnData] = None,
        to: str = "collectri",
        organism="human",
        base_pr_threshold=0,
    ):
        if other is None:
            print("loading GT, ", to)
            gt = get_GT_db(name=to, organism=organism)
        else:
            other = other.grn[(other.grn.sum(0) != 0)]
            other = other.loc[:, other.sum() > 0]
            gt = []
            for k, v in zip(*np.where(other != 0)):
                gt.append([other.index[k], other.columns[v], 1])
            gt = pd.DataFrame(gt, columns=["source", "target", "weight"])

        varnames = list(set(gt.iloc[:, :2].values.flatten()))
        intersection = self.grn.grn.index.intersection(varnames)

        gt = gt[gt.source.isin(intersection) & gt.target.isin(intersection)]
        true_con = set([frozenset([i, j]) for i, j in gt.iloc[:, :2].values])

        varnames = list(set(gt.iloc[:, :2].values.flatten()))
        if self.doplot:
            print("intersection of {} genes".format(intersection.size))
            print("intersection pct:", intersection.size / self.grn.grn.index.size)
        grn = self.grn.grn.loc[varnames, varnames]
        return compute_auprc(
            grn, true_con, base_pr_threshold=base_pr_threshold, doplot=self.doplot
        )

    def scprint_benchmark(self, base_pr_threshold=0):
        print("base enrichment")
        metrics = {}
        for elem in ["Central", "Targets"]:
            res = utils.enrichment(
                self.grn,
                of=elem,
                gene_sets=[
                    {"TFs": utils.TF},
                    utils.file_dir + "/celltype.gmt",
                ],
                doplot=False,
                top_k=10,
            )
            if (
                len(res.res2d[(res.res2d["FDR q-val"] < 0.1) & (res.res2d["NES"] > 1)])
                > 0
            ):
                metrics.update(
                    {
                        "enriched_terms_"
                        + elem: res.res2d[
                            (res.res2d["FDR q-val"] < 0.1) & (res.res2d["NES"] > 1)
                        ].Term.tolist()
                    }
                )
                if self.doplot:
                    _ = res.plot(
                        terms=res.res2d[
                            (res.res2d["FDR q-val"] < 0.1) & (res.res2d["NES"] > 1)
                        ]
                        .sort_values(by=["NES"], ascending=False)
                        .Term.iloc[0]
                    )
            if self.doplot:
                try:
                    _ = res.plot(terms="0__TFs")
                    plt.show()
                except KeyError:
                    pass
        print("_________________________________________")
        print("TF specific enrichment")
        tfchip = gp.get_library(name="ENCODE_TF_ChIP-seq_2014")
        TFinchip = {i: i.split(" ")[0] for i in tfchip.keys()}
        res = {}
        i, j = 0, 0
        for k, v in TFinchip.items():
            if v not in self.grn.grn.columns:
                continue
            # print(k)
            j += 1
            test = self.grn.grn.loc[v].sort_values(ascending=False)
            if len(set(test.index) & set(tfchip[k])) == 0:
                continue
            if test.sum() == 0:
                continue
            pre_res = gp.prerank(
                rnk=test,  # or rnk = rnk,
                gene_sets=[
                    # "Chromosome_Location",
                    {v: tfchip[k]}
                ],
                min_size=1,
                max_size=4000,
                permutation_num=1000,
            )
            val = (
                pre_res.res2d[
                    (pre_res.res2d["FDR q-val"] < 0.1) & (pre_res.res2d["NES"] > 1)
                ]
                .sort_values(by=["NES"], ascending=False)
                .drop(columns=["Name"])
            )
            if len(val.Term.tolist()) > 0:
                print("found! ", val.Term.tolist()[0])
                print(pre_res.res2d["NES"])
                print("\n")
                i += 1
            res[k] = pre_res.res2d
            # print(val.Term.tolist()[:2])
        print("found some significant results for ", i * 100 / j, "% TFs\n")
        metrics.update({"significant_enriched_TFtargets": i * 100 / j})
        print("_________________________________________")
        metrics.update(self.compare_to(to="omnipath", organism="human"))
        return metrics


def train_classifier(
    adj,
    genes,
    gt="omnipath",
    other=None,
    train_size=1_000_000,
    doplot=True,
    class_weight={1: 200, 0: 1},
):
    if other is not None:
        pass
    else:
        gt = get_GT_db(name=gt)

    varnames = set(gt.iloc[:, :2].values.flatten())
    intersection = varnames & set(genes)

    loc = np.isin(genes, np.array(list(intersection)))
    genes = np.array(genes)[loc].tolist()
    adj = adj[:, loc, :][loc, :, :]

    da = np.zeros((len(genes), len(genes)), dtype=np.float)
    for i, j in gt.iloc[:, :2].values:
        if i in genes and j in genes:
            da[genes.index(i), genes.index(j)] = 1
            da[genes.index(j), genes.index(i)] = 1
    print("true elem", da.sum())
    da = da.flatten()
    adj = adj.reshape(-1, adj.shape[-1])

    X_train, X_test, y_train, y_test = train_test_split(
        adj, da, random_state=0, train_size=train_size
    )

    clf = LogisticRegression(
        penalty="l1",
        C=1,
        solver="liblinear",
        class_weight=class_weight,
        max_iter=10_000,
        n_jobs=8,
    )
    clf.fit(X_train, y_train)
    print("coef", clf.coef_)

    pred = clf.predict(X_test)
    print("precision", (pred[y_test == 1] == 1).sum() / (pred == 1).sum())
    print("random precision", y_test.sum() / len(y_test))
    print("recall", (pred[y_test == 1] == 1).sum() / y_test.sum())
    print("predicted true", pred.sum())
    print("number of true", y_test.sum())
    if doplot:
        PrecisionRecallDisplay.from_estimator(
            clf, X_test, y_test, plot_chance_level=True
        )
        plt.show()


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
                + "/../data/GroundTruth/stone_and_sroy/scRNA/liu_rna_filtered_log2.tsv.gz",
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
                + "/../data/GroundTruth/stone_and_sroy/scRNA/chen_rna_filtered_log2.tsv.gz",
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
    adata = concat([adata_liu, adata_chen], join=join, fill_value=0)
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
    grn.var["TFs"] = [True if i in utils.TF else False for i in grn.var_names]
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
    grn = GRNAnnData(grn=VIM, X=mat, var=var, obs=adata.obs)
    grn.var_names = grn.var["symbol"]
    grn.var["TFs"] = [True if i in utils.TF else False for i in grn.var_names]
    return grn


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
    elif name == "omnipath":
        if not os.path.exists(FILEDIR + "/../data/omnipath.parquet"):
            interactions = AllInteractions()
            net = interactions.get(exclude=["small_molecule", "lncrna_mrna"])
            hgnc = Annotations.get(resources="HGNC")
            rename = {v.uniprot: v.genesymbol for k, v in hgnc.iterrows()}
            net.source = net.source.replace(rename)
            net.target = net.target.replace(rename)
            net.to_parquet(FILEDIR + "/../data/omnipath.parquet")
        else:
            net = pd.read_parquet(FILEDIR + "/../data/omnipath.parquet")
    else:
        raise ValueError(f"provided name: '{name}' is not amongst the available names.")
    # varnames = list(set(net.iloc[:, :2].values.flatten()))
    # adata = AnnData(var=varnames)
    # adata.var_names = varnames
    # grn = from_adata_and_longform(adata, net, has_weight=True)
    # return grn
    return net


def pd_load_cached(url, loc="/tmp/", cache=True, **kwargs):
    # Check if the file exists, if not, download it
    loc += url.split("/")[-1]
    if not os.path.isfile(loc) or not cache:
        urllib.request.urlretrieve(url, loc)
    # Load the data from the file
    return pd.read_csv(loc, **kwargs)


def compute_auprc(grn: pd.DataFrame, true_con: set, base_pr_threshold=0, doplot=True):
    metrics = {}
    grn_con = set(
        [
            frozenset([grn.index[i], grn.columns[j]])
            for i, j in zip(*np.where(grn > base_pr_threshold))
        ]
    )
    rand = grn.shape[0] * grn.shape[1]  # / 2) + (grn.shape[0] / 2))
    precision = len(true_con & grn_con) / len(grn_con)
    recall = len(true_con & grn_con) / len(true_con)
    rand_rec = len(grn_con) / rand
    rand_prec = len(true_con) / rand

    if doplot:
        print(
            "precision: ",
            precision,
            "\nrecall: ",
            recall,
            "\nrandom recall:",
            rand_rec,
            "\nrandom precision:",
            rand_prec,
        )
    # Initialize lists to store precision and recall values
    precision_list = [precision]
    recall_list = [recall]
    rand_recall_list = [rand_rec]
    # Define the thresholds to vary
    thresholds = np.linspace(grn.values.min(), grn.values.max(), 50)
    # Calculate precision and recall for each threshold
    for threshold in tqdm.tqdm(thresholds[1:]):
        grn_con_threshold = set(
            [
                frozenset([grn.index[i], grn.columns[j]])
                for i, j in zip(*np.where(grn >= threshold))
            ]
        )
        precision = len(true_con & grn_con_threshold) / len(grn_con_threshold)
        recall = len(true_con & grn_con_threshold) / len(true_con)
        rand_rec = len(grn_con_threshold) / rand
        precision_list.append(precision)
        recall_list.append(recall)
        rand_recall_list.append(rand_rec)

    # Calculate AUPRC by integrating the precision-recall curve
    auprc = -np.trapz(precision_list, recall_list)
    metrics["auprc"] = auprc
    maxv = np.max(np.array(recall_list) - np.array(rand_recall_list))
    meanv = np.mean(np.array(recall_list) - np.array(rand_recall_list))
    metrics.update({"pr_increase_to_random": (meanv, maxv)})
    if doplot:
        print("Area Under Precision-Recall Curve (AUPRC): ", auprc)
        print("overal increase: (mean, max)", (meanv, maxv))
        # compute EPR
        # Flatten the array and get the indices of the top K values
        indices = np.argpartition(grn.values.flatten(), -len(true_con))[
            -len(true_con) :
        ]
        # Convert the indices to 2D
        top_K_indices = np.column_stack(np.unravel_index(indices, grn.shape))
        grn_con_topk = set(
            [frozenset([grn.index[i], grn.columns[j]]) for i, j in top_K_indices]
        )
        # Compute the odds ratio
        true_positive = len(true_con & grn_con_topk)
        false_positive = len(grn_con_topk) - true_positive
        false_negative = len(true_con) - true_positive
        true_negative = rand - true_positive - false_positive - false_negative
        print("true_positive", true_positive)
        print("True Negative: ", true_negative)
        print("False Positive: ", false_positive)
        print("False Negative: ", false_negative)

        # Avoid division by zero
        if true_negative == 0 or false_positive == 0:
            odds_ratio = float("inf")
        else:
            odds_ratio = (true_positive * true_negative) / (
                false_positive * false_negative
            )

        print("Odds Ratio: ", odds_ratio)

        plt.figure(figsize=(10, 8))
        plt.plot(
            recall_list,
            precision_list,
            marker=".",
            linestyle="-",
            color="b",
            label="p-r",
        )
        plt.plot(
            rand_recall_list,
            precision_list,
            marker=".",
            linestyle="-",
            color="r",
            label="rand rec",
        )
        plt.legend(loc="lower left")
        plt.title("Precision-Recall Curve")
        plt.xlabel("Recall")
        plt.ylabel("Precision")
        plt.xscale("log")
        plt.grid(True)
        plt.show()
    return metrics
