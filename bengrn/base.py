"""
bengrn base module.
"""

import pandas as pd
import urllib.request
import os.path
import json

from grnndata import GRNAnnData, from_adata_and_longform, from_scope_loomfile, utils
from anndata import AnnData, concat
from typing import Optional, Union
from .tools import GENIE3
import numpy as np
from anndata.utils import make_index_unique

from arboreto.algo import grnboost2

# issue here of using an older version of numpy calling np.object instead of np.object_
# from pyscenic.utils import modules_from_adjacencies
# from pyscenic.prune import prune2df, df2regulons
# from pyscenic.aucell import aucell

from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase

import scipy.stats
from sklearn.metrics import precision_recall_curve, auc, PrecisionRecallDisplay
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.linear_model import RidgeClassifier

import gseapy as gp

import tqdm
import logging

import matplotlib.pyplot as plt

from scipy.sparse import csr_matrix, csc_matrix

import scanpy as sc

# example constant variable
NAME = "bengrn"

FILEDIR = os.path.dirname(os.path.realpath(__file__))


# Define a function to calculate precision and recall
def precision_recall(true_con, grn_con):
    """
    Calculate precision and recall from the true and predicted connections.

    Args:
        true_con (np.array): The true connections.
        grn_con (np.array): The predicted connections.

    Returns:
        float: The precision value.
        float: The recall value.
    """
    tp = len(true_con & grn_con)
    fp = len(grn_con - true_con)
    fn = len(true_con - grn_con)

    precision = tp / (tp + fp)
    recall = tp / (tp + fn)

    return precision, recall


class BenGRN:
    def __init__(
        self,
        grn: GRNAnnData,
        full_dataset: Optional[AnnData] = None,
        doplot: bool = True,
        do_auc: bool = True,
    ):
        """
        Initializes the BenGRN class.

        Args:
            grn (GRNAnnData): The Gene Regulatory Network data.
            full_dataset (Optional[AnnData]): The full dataset, defaults to None.
            doplot (bool): Whether to plot the results, defaults to True.
            do_auc (bool): Whether to calculate the Area Under the Precision-Recall Curve, defaults to True.
        """
        self.grn = grn
        self.full_dataset = full_dataset
        self.doplot = doplot
        self.do_auc = do_auc

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
        organism: str = "human",
    ):
        """
        compare_to compares the GRN to another GRN.

        Args:
            other (Optional[GRNAnnData], optional): The other GRN to compare to. Defaults to None.
                If not given can use a default GRN from the 'to' argument.
            to (str, optional): The name of the other GRN to compare to. Defaults to "collectri".
                If 'other' is given, this argument is ignored.
            organism (str, optional): The organism of the GRN to compare to. Defaults to "human".

        Returns:
            dict: The metrics of the comparison.
        """
        if other is None:
            if self.doplot:
                print("loading GT, ", to)
            gt = get_GT_db(name=to, organism=organism)
            # gt = gt[gt.type != "post_translational"]
            varnames = set(gt.iloc[:, :2].values.flatten())
            intersection = varnames & set(self.grn.var["symbol"].tolist())
            loc = self.grn.var["symbol"].isin(intersection)
            adj = self.grn.varp["GRN"][:, loc][loc, :]
            genes = self.grn.var.loc[loc, "symbol"].tolist()

            da = np.zeros(adj.shape, dtype=float)
            for i, j in gt.iloc[:, :2].values:
                if i in genes and j in genes:
                    da[genes.index(i), genes.index(j)] = 1
            if self.doplot:
                print("intersection of {} genes".format(len(intersection)))
                print("intersection pct:", len(intersection) / len(self.grn.grn.index))
        else:
            elems = other.var[other.grn.sum(1) != 0].index.tolist()
            da = other.get(self.grn.var.index.tolist()).get(elems).targets
            if da.shape[1] < 5:
                print("da is very small: ", da.shape[1])
            # da = da.iloc[6:]
            adj = self.grn.grn.loc[da.index.values, da.columns.values].values
            da = da.values

        return compute_pr(
            adj,
            da,
            doplot=self.doplot,
            do_auc=self.do_auc,
        )

    def scprint_benchmark(self, elems=["Central", "Regulators", "Targets"]):
        """
        scprint_benchmark full benchmarks of the GRN as in the scPRINT paper.

        It will apply first an enrichment analysis over the [elems] of the GRN looking for TF enrichment and cell type marker gene enrichment
        It will then apply an enrichment over each TF in the GRN for their targets in ENCODE.
        Finaly, it will compare it to the OmniPath database GRN using precision recall type metrics.

        Args:
            elems (list, optional): The genes in the GRN, to benchmark. Defaults to ["Central", "Regulators", "Targets"].
                It corresponds to different views of the network.

        Returns:
            dict: The metrics of the benchmark.
        """
        print("base enrichment")
        metrics = {}
        for elem in elems:
            if elem == "Central" and (self.grn.varp["GRN"] != 0).sum() > 100_000_000:
                print("too many genes for central computation")
                continue
            res = utils.enrichment(
                self.grn,
                of=elem,
                gene_sets=[
                    {"TFs": [i.split(".")[0] for i in utils.TF]},
                    utils.file_dir + "/celltype.gmt",
                ],
                doplot=False,
                maxsize=2000,
                top_k=10,
            )
            if res is None:
                continue
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
                try:
                    _ = res.plot(terms="0__TFs")
                    plt.show()
                except KeyError:
                    pass
            istrue = metrics.get("TF_enr", False)
            istrue = istrue or (
                res.res2d.loc[res.res2d.Term == "0__TFs", "FDR q-val"].iloc[0] < 0.1
            )
            metrics.update({"TF_enr": istrue})
        if self.doplot:
            print("_________________________________________")
            print("TF specific enrichment")
        with open(FILEDIR + "/../data/tfchip_data.json", "r") as file:
            tfchip = json.load(file)
        TFinchip = {i: i.split(" ")[0] for i in tfchip.keys()}
        res = {}
        i, j = 0, 0
        previous_level = logging.root.manager.disable
        logging.disable(logging.WARNING)
        for k, v in TFinchip.items():
            if v not in self.grn.grn.columns:
                continue
            # print(k)
            j += 1
            test = self.grn.grn.T.loc[[v]].sort_values(by=v, ascending=False).T
            if len(set(test.index) & set(tfchip[k])) == 0:
                continue
            if test.iloc[:, 0].sum() == 0:
                continue
            try:
                pre_res = gp.prerank(
                    rnk=test,
                    gene_sets=[{v: tfchip[k]}],
                    background=self.grn.var.index.tolist(),
                    min_size=1,
                    max_size=4000,
                    permutation_num=1000,
                )
            except IndexError:
                continue
            val = (
                pre_res.res2d[
                    (pre_res.res2d["FDR q-val"] < 0.05) & (pre_res.res2d["NES"] > 1)
                ]
                .sort_values(by=["NES"], ascending=False)
                .drop(columns=["Name"])
            )
            if len(val.Term.tolist()) > 0:
                # print("found! ", val.Term.tolist()[0])
                # print(pre_res.res2d["NES"])
                # print(pre_res.res2d["FDR q-val"])
                # print("\n")
                i += 1
            else:
                pass
                # print("no sig...")
            res[k] = pre_res.res2d
            # print(val.Term.tolist()[:2])
        logging.disable(previous_level)
        j = j if j != 0 else 1
        if self.doplot:
            print("found some significant results for ", i * 100 / j, "% TFs\n")
            print("_________________________________________")
        metrics.update({"significant_enriched_TFtargets": i * 100 / j})
        metrics.update(self.compare_to(to="omnipath", organism="human"))
        return metrics


def train_classifier(
    grn: GRNAnnData,
    gt: str = "omnipath",
    other: Optional[GRNAnnData] = None,
    use_col: str = "symbol",
    train_size: float = 0.2,
    doplot: bool = True,
    class_weight: dict = {1: 200, 0: 1},
    max_iter: int = 1_000,
    C: float = 1.0,
    return_full: bool = True,
    shuffle: bool = False,
    **kwargs,
):
    """
    train_classifier trains a classifier to generate a GRN that maps to the ground truth.

    Uses a RidgeClassifier to select the best combination of networks to predict the ground truth.
    It is used for the head classification part in the scPRINT paper.

    Args:
        grn (GRNAnnData): The Gene Regulatory Network data.
        gt (str, optional): The name of the ground truth database to use. Defaults to "omnipath".
        other (GRNAnnData, optional): Another GRN to compare against. Defaults to None.
        use_col (str, optional): The column name to use for gene symbols. Defaults to "symbol".
        train_size (float, optional): The proportion of the dataset to include in the train split. Defaults to 0.2.
        doplot (bool, optional): Whether to plot the results. Defaults to True.
        class_weight (dict, optional): Weights associated with classes in the form {class_label: weight}. Defaults to {1: 200, 0: 1}.
        max_iter (int, optional): Maximum number of iterations for the classifier. Defaults to 1_000.
        C (float, optional): Regularization strength; must be a positive float. Defaults to 1.0.
        return_full (bool, optional): Whether to return the full classifier object. Defaults to True.
        shuffle (bool, optional): Whether to shuffle the data before splitting. Defaults to False.

    Returns:
        (GRNAnnData, dict, RidgeClassifier): The Gene Regulatory Network data, the metrics of the classifier, and the classifier object.
    """
    if other is not None:
        elems = other.var[other.grn.sum(1) != 0].index.tolist()
        sub = other.get(grn.var[use_col].tolist()).get(elems).targets
        if sub.shape[1] < 5:
            print("sub is very small: ", sub.shape[1])
        genes = grn.var[use_col].tolist()
        args = np.argsort(genes)
        genes = np.array(genes)[args]
        adj = grn.varp["GRN"][args, :, :][:, args, :][np.isin(genes, sub.index.values)][
            :, np.isin(genes, sub.columns.values)
        ]
        print("pred shape", adj.shape)
        da = sub.values
    else:
        gt = get_GT_db(name=gt)
        varnames = set(gt.iloc[:, :2].values.flatten())
        intersection = varnames & set(grn.var[use_col].tolist())
        loc = grn.var[use_col].isin(intersection)
        adj = grn.varp["GRN"][:, loc, :][loc, :, :]
        genes = grn.var.loc[loc][use_col].tolist()

        da = np.zeros((len(genes), len(genes)), dtype=float)
        for i, j in gt.iloc[:, :2].values:
            if i in genes and j in genes:
                da[genes.index(i), genes.index(j)] = 1

    print("true elem", int(da.sum()), "...")
    da = da.flatten()
    adj = adj.reshape(-1, adj.shape[-1])

    X_train, X_test, y_train, y_test = train_test_split(
        adj, da, random_state=0, train_size=train_size, shuffle=shuffle
    )
    print("doing classification....")

    clf = RidgeClassifier(
        alpha=C,
        fit_intercept=False,
        class_weight=class_weight,
        # solver="saga",
        max_iter=max_iter,
        positive=True,
    )
    # clf = LogisticRegression(
    #    penalty="l1",
    #    C=C,
    #    solver="saga",
    #    class_weight=class_weight,
    #    max_iter=max_iter,
    #    n_jobs=8,
    #    fit_intercept=False,
    #    verbose=10,
    #    **kwargs,
    # )
    clf.fit(X_train, y_train)
    pred = clf.predict(X_test)
    # epr = compute_epr(clf, X_test, y_test)
    metrics = {
        "used_heads": (clf.coef_ != 0).sum(),
        "precision": (pred[y_test == 1] == 1).sum() / (pred == 1).sum(),
        "random_precision": y_test.sum() / len(y_test),
        "recall": (pred[y_test == 1] == 1).sum() / y_test.sum(),
        "predicted_true": pred.sum(),
        "number_of_true": y_test.sum(),
        # "epr": epr,
    }
    if doplot:
        print("metrics", metrics)
        PrecisionRecallDisplay.from_estimator(
            clf, X_test, y_test, plot_chance_level=True
        )
        plt.show()
    if return_full:
        adj = grn.varp["GRN"]
        grn.varp["classified"] = clf.predict(adj.reshape(-1, adj.shape[-1])).reshape(
            len(grn.var), len(grn.var)
        )
    return grn, metrics, clf


def get_scenicplus(
    filepath=FILEDIR + "/../data/10xPBMC_homo_scenicplus_genebased_scope.loom",
):
    """
    This function retrieves a loomx scenicplus data from a given file path and loads it as a GrnnData

    Args:
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


def get_sroy_gt(
    get: str = "main", join: str = "outer", species: str = "human", gt: str = "full"
) -> GRNAnnData:
    """
    This function retrieves the ground truth data from the McCall et al.'s paper.

    Args:
        get (str): The specific dataset to retrieve. Options include "main", "liu", and "chen".
        join (str, optional): The type of join to be performed when concatenating the data. Default is "outer".
        species (str, optional): The species of the dataset. Default is "human".
        gt (str, optional): The type of ground truth data to retrieve. Options include "full", "chip", and "ko". Default is "full".

    Returns:
        GrnAnnData: The ground truth data as a grnndata object
    """
    if species == "human":
        if gt == "full":
            df = pd.read_csv(
                FILEDIR + "/../data/GroundTruth/stone_and_sroy/hESC_ground_truth.tsv",
                sep="\t",
                header=None,
            )
        elif gt == "chip":
            df = pd.read_csv(
                FILEDIR
                + "/../data/GroundTruth/stone_and_sroy/gold_standards/hESC/hESC_chipunion.txt",
                sep="\t",
                header=None,
            )
        elif gt == "ko":
            df = pd.read_csv(
                FILEDIR
                + "/../data/GroundTruth/stone_and_sroy/gold_standards/hESC/hESC_KDUnion.txt",
                sep="\t",
                header=None,
            )
        if get == "liu":
            adata = AnnData(
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
        elif get == "chen":
            adata = AnnData(
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
        elif get == "han":
            adata = AnnData(
                unnormalize(
                    pd.read_csv(
                        FILEDIR
                        + "/../data/GroundTruth/stone_and_sroy/scRNA/human_han_GSE107552.csv.gz",
                    )
                    .set_index("Cell")
                    .T,
                    is_root=True,
                )
            )
        elif get == "mine":
            # https://www.ebi.ac.uk/gxa/sc/experiments/E-GEOD-36552/downloads
            adata = sc.read_mtx(
                FILEDIR
                + "/../data/GroundTruth/stone_and_sroy/scRNA/E-GEOD-36552.aggregated_filtered_counts.mtx"
            ).T
            col = pd.read_csv(
                FILEDIR
                + "/../data/GroundTruth/stone_and_sroy/scRNA/E-GEOD-36552.aggregated_filtered_counts.mtx_rows",
                header=None,
                sep="\t",
            )
            adata.var.index = col[0]
            genesdf = load_genes()
            intersect_genes = set(adata.var.index).intersection(set(genesdf.index))
            adata = adata[:, list(intersect_genes)]
            genesdf = genesdf.loc[adata.var.index]
            adata.var["ensembl_id"] = adata.var.index
            adata.var.index = make_index_unique(genesdf["symbol"].astype(str))
        else:
            raise ValueError("get must be one of 'liu', 'chen', 'han', or 'mine'")
        adata.obs["organism_ontology_term_id"] = "NCBITaxon:9606"
    elif species == "mouse":
        if gt == "full":
            df = pd.read_csv(
                FILEDIR + "/../data/GroundTruth/stone_and_sroy/mESC_ground_truth.tsv",
                sep="\t",
                header=None,
            )
        elif gt == "chip":
            df = pd.read_csv(
                FILEDIR
                + "/../data/GroundTruth/stone_and_sroy/gold_standards/mESC/mESC_chipunion.txt",
                sep="\t",
                header=None,
            )
        elif gt == "ko":
            df = pd.read_csv(
                FILEDIR
                + "/../data/GroundTruth/stone_and_sroy/gold_standards/mESC/mESC_KDUnion.txt",
                sep="\t",
                header=None,
            )
        if get == "duren":
            adata = AnnData(
                (
                    2
                    ** pd.read_csv(
                        FILEDIR
                        + "/../data/GroundTruth/stone_and_sroy/scRNA/duren_rna_filtered_log2.tsv.gz",
                        sep="\t",
                    )
                )
                - 1
            ).T
        elif get == "semrau":
            adata = AnnData(
                (
                    2
                    ** pd.read_csv(
                        FILEDIR
                        + "/../data/GroundTruth/stone_and_sroy/scRNA/semrau_rna_filtered_log2.tsv.gz",
                        sep="\t",
                    )
                )
                - 1
            ).T
        elif get == "tran":
            adata = AnnData(
                unnormalize(
                    pd.read_csv(
                        FILEDIR
                        + "/../data/GroundTruth/stone_and_sroy/scRNA/mouse_tran_A2S.csv.gz",
                    )
                    .set_index("Cell")
                    .T,
                    is_root=True,
                )
            )
        elif get == "zhao":
            adata = AnnData(
                unnormalize(
                    pd.read_csv(
                        FILEDIR
                        + "/../data/GroundTruth/stone_and_sroy/scRNA/mouse_zhao_GSE114952.csv.gz",
                    )
                    .set_index("Cell")
                    .T,
                    is_root=True,
                )
            )
        else:
            raise ValueError("get must be one of 'duren', 'semrau', 'tran', or 'zhao'")
        adata.obs["organism_ontology_term_id"] = "NCBITaxon:10090"
    return from_adata_and_longform(adata, df)


def get_perturb_gt(
    url_bh: str = "https://plus.figshare.com/ndownloader/files/38349308",
    filename_bh: str = FILEDIR + "/../data/BH-corrected.csv.gz",
    url_adata: str = "https://plus.figshare.com/ndownloader/files/35773219",
    filename_adata: str = FILEDIR + "/../data/ess_perturb_sc.h5ad",
):
    """
    get_perturb_gt retrieves the genome wide perturb seq ground truth data.

    Args:
        url_bh (str, optional): The URL to download the BH-corrected data. Defaults to "https://plus.figshare.com/ndownloader/files/38349308".
        filename_bh (str, optional): The local filename to save the BH-corrected data. Defaults to FILEDIR + "/../data/BH-corrected.csv.gz".
        url_adata (str, optional): The URL to download the single-cell perturbation data. Defaults to "https://plus.figshare.com/ndownloader/files/35773219".
        filename_adata (str, optional): The local filename to save the single-cell perturbation data. Defaults to FILEDIR + "/../data/ess_perturb_sc.h5ad".

    Returns:
        GRNAnnData: The Gene Regulatory Network data as a GRNAnnData object.
    """
    if not os.path.exists(filename_bh):
        urllib.request.urlretrieve(url_bh, filename_bh)
    pert = pd.read_csv(filename_bh)
    pert = pert.set_index("Unnamed: 0").T
    pert.index = [i.split("_")[-1] for i in pert.index]
    pert = pert[~pert.index.duplicated(keep="first")].T
    pert = pert < 0.05
    adata_sc = sc.read(
        filename_adata,
        backup_url=url_adata,
    )
    adata_sc = adata_sc[adata_sc.obs.gene_id == "non-targeting"]
    adata_sc[:, adata_sc.var.index.isin(set(pert.index) | set(pert.columns))]
    adata_sc.obs["organism_ontology_term_id"] = "NCBITaxon:9606"
    adata_sc = adata_sc[:, adata_sc.var.sort_index().index]

    pert = pert.loc[pert.index.isin(adata_sc.var.index)].loc[
        :, pert.columns.isin(adata_sc.var.index)
    ]

    missing_indices = list(set(adata_sc.var.index) - set(pert.index))
    pert = pert.reindex(pert.index.union(missing_indices), fill_value=False)
    missing_indices = list(set(adata_sc.var.index) - set(pert.columns))
    pert = pert.reindex(columns=pert.columns.union(missing_indices), fill_value=False)

    pert = pert.loc[adata_sc.var.index].loc[:, adata_sc.var.index]
    return GRNAnnData(
        X=csr_matrix(adata_sc.X.toarray()),
        var=adata_sc.var,
        obs=adata_sc.obs,
        grn=csr_matrix(pert.values),
    )


# def compute_scenic(adata, data_dir=FILEDIR + "/../data"):
#    """
#    This function computes the SCENIC algorithm on the given dataset.
#
#    Args:
#        adata (AnnData): The annotated data matrix of shape n_obs x n_vars. Rows correspond to cells and columns to genes.
#        data_dir (str, optional): The directory where the data files will be stored. Defaults to FILEDIR + "/../data".
#
#    Returns:
#        GRNAnnData: The Gene Regulatory Network data.
#    """
#    raise NotImplementedError("This function is not implemented yet")
#    os.makedirs(data_dir, exist_ok=True)
#
#    url1 = "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
#    url2 = "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
#
#    file1 = os.path.join(data_dir, os.path.basename(url1))
#    file2 = os.path.join(data_dir, os.path.basename(url2))
#
#    if not os.path.exists(file1):
#        urllib.request.urlretrieve(url1, file1)
#    if not os.path.exists(file2):
#        urllib.request.urlretrieve(url2, file2)
#    dbs = RankingDatabase(fname=file2, name="genes_10kb_human_rank_v10")
#    adjacencies = grnboost2(
#        expression_data=adata.to_df(), tf_names=utils.TF, verbose=True
#    )
#    modules = list(modules_from_adjacencies(adjacencies, adata.to_df()))
#
#    df = prune2df([dbs], modules, file1)
#    regulons = df2regulons(df)
#    # auc_mtx = aucell(adata.to_df(), regulons, num_workers=8)
#    # sns.clustermap(auc_mtx, figsize=(12,12))
#
#    # compute the grn matrix
#    var_names = adata.var_names.tolist()
#    da = np.zeros((len(var_names), len(var_names)), dtype=float)
#    for reg in regulons:
#        i = var_names.index(reg.transcription_factor)
#        for k, v in reg.gene2weight.items():
#            da[i, var_names.index(k)] = v
#    # create the grndata
#    grn = GRNAnnData(adata, grn=da)
#
#    da = np.zeros((len(var_names), len(var_names)), dtype=float)
#    for reg in modules:
#        i = var_names.index(reg.transcription_factor)
#        for k, v in reg.gene2weight.items():
#            da[i, var_names.index(k)] = v
#    grn.varp["modules"] = da
#    grn.var["TFs"] = [True if i in utils.TF else False for i in grn.var_names]
#    return grn


def compute_genie3(
    adata: AnnData, nthreads: int = 30, ntrees: int = 100, **kwargs
) -> GRNAnnData:
    """
    This function computes the GENIE3 algorithm on the given data.

    Args:
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


def get_GT_db(
    name: str = "collectri", organism: str = "human", split_complexes: bool = True
) -> pd.DataFrame:
    """
    use_prior_network loads a prior GRN from a list of available networks.

    Args:
        name (str, optional): name of the network to load. Defaults to "collectri".
        organism (str, optional): organism to load the network for. Defaults to "human".
        split_complexes (bool, optional): whether to split complexes into individual genes. Defaults to True.

    Returns:
        pd.DataFrame: The prior GRN as a pandas DataFrame.

    Raises:
        ValueError: if the provided name is not amongst the available names.
    """
    if name == "tflink":
        TFLINK = "https://cdn.netbiol.org/tflink/download_files/TFLink_Homo_sapiens_interactions_All_simpleFormat_v1.0.tsv.gz"
        net = pd_load_cached(TFLINK)
        net = net.rename(columns={"Name.TF": "source", "Name.Target": "target"})
    elif name == "htftarget":
        HTFTARGET = "http://bioinfo.life.hust.edu.cn/static/hTFtarget/file_download/tf-target-infomation.txt"
        net = pd_load_cached(HTFTARGET)
        net = net.rename(columns={"TF": "source"})
    elif name == "collectri":
        import decoupler as dc

        net = dc.get_collectri(organism=organism, split_complexes=split_complexes).drop(
            columns=["PMID"]
        )
    elif name == "dorothea":
        import decoupler as dc

        net = dc.get_dorothea(organism=organism)
    elif name == "omnipath":
        if not os.path.exists(FILEDIR + "/../data/omnipath.parquet"):
            from omnipath.interactions import AllInteractions
            from omnipath.requests import Annotations

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


def compute_pr(
    grn: np.array,
    true: np.array,
    do_auc: bool = True,
    doplot: bool = True,
):
    """
    compute_pr computes the precision and recall metrics for the given GRN and true matrix.

    Args:
        grn (np.array): The Gene Regulatory Network matrix, where each element represents the strength of the regulatory relationship between genes.
        true (np.array): The ground truth matrix, where each element indicates the presence (1) or absence (0) of a regulatory relationship.
        do_auc (bool, optional): Whether to compute the Area Under the Precision-Recall Curve (AUPRC). Defaults to True.
        doplot (bool, optional): Whether to plot the precision and recall metrics. Defaults to True.

    Raises:
        ValueError: If the shape of the GRN and the true matrix do not match.

    Returns:
        dict: A dictionary containing precision, recall, and random precision metrics.
    """
    if grn.shape != true.shape:
        raise ValueError("The shape of the GRN and the true matrix do not match.")
    metrics = {}
    if isinstance(grn, (csr_matrix, csc_matrix)):
        grn = grn.toarray()
    if isinstance(true, (csr_matrix, csc_matrix)):
        true = true.toarray()
    true = true.astype(bool)
    tot = (grn.shape[0] * grn.shape[1]) - grn.shape[0]
    precision = (grn[true] != 0).sum() / (grn != 0).sum()
    recall = (grn[true] != 0).sum() / true.sum()
    rand_prec = true.sum() / tot

    if doplot:
        print(
            "precision: ",
            precision,
            "\nrecall: ",
            recall,
            "\nrandom precision:",
            rand_prec,
        )
    metrics.update(
        {
            "precision": precision,
            "recall": recall,
            "rand_precision": rand_prec,
        }
    )
    # Initialize lists to store precision and recall values
    precision_list = [precision]
    recall_list = [recall]
    # Define the thresholds to vary
    thresholds = np.append(
        np.linspace(0, 1, 101)[:-2], np.log10(np.logspace(0.99, 1, 30))
    )
    thresholds = np.quantile(grn, thresholds)
    # Calculate precision and recall for each threshold
    if do_auc:
        for threshold in tqdm.tqdm(thresholds[1:]):
            precision = (grn[true] > threshold).sum() / (grn > threshold).sum()
            recall = (grn[true] > threshold).sum() / true.sum()
            precision_list.append(precision)
            recall_list.append(recall)
        # Calculate AUPRC by integrating the precision-recall curve
        if 1.0 not in recall_list:
            precision_list.insert(0, rand_prec)
            recall_list.insert(0, recall_list[0])
            precision_list.insert(0, rand_prec)
            recall_list.insert(0, 1.0)
        precision_list = np.nan_to_num(np.array(precision_list))
        recall_list = np.nan_to_num(np.array(recall_list))
        auprc = -np.trapz(precision_list, recall_list)
        metrics["auprc"] = auprc

        # Compute Average Precision (AP) manually
        sorted_indices = np.argsort(-grn.flatten())
        sorted_true = true.flatten()[sorted_indices]

        tp_cumsum = np.cumsum(sorted_true)
        fp_cumsum = np.cumsum(~sorted_true)

        precision_at_k = tp_cumsum / (tp_cumsum + fp_cumsum)
        recall_at_k = tp_cumsum / true.sum()

        ap = np.sum(precision_at_k[1:] * np.diff(recall_at_k))
        metrics["ap"] = ap
        if doplot:
            print("Average Precision (AP): ", ap)
        if doplot:
            print("Area Under Precision-Recall Curve (AUPRC): ", auprc)

    # compute EPR
    # get the indices of the topK highest values in "grn"
    if isinstance(grn, csr_matrix):
        grn = grn.toarray()
    if isinstance(grn, csc_matrix):
        grn = grn.toarray()
    indices = np.argpartition(grn.flatten(), -int(true.sum()))[-int(true.sum()) :]
    # Compute the odds ratio
    true_positive = true[np.unravel_index(indices, true.shape)].sum()
    false_positive = true.sum() - true_positive
    # this is normal as we compute on the same number of pred_pos as true_pos
    false_negative = true.sum() - true_positive
    true_negative = tot - true_positive - false_positive - false_negative
    # Avoid division by zero
    # this is a debugger line
    if true_negative == 0 or false_positive == 0:
        odds_ratio = float("inf")
    else:
        odds_ratio = (true_positive * true_negative) / (false_positive * false_negative)

    metrics.update({"epr": odds_ratio})
    if doplot:
        print("EPR:", odds_ratio)
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
            [recall_list[0], recall_list[-1]],
            [rand_prec, rand_prec],
            linestyle="--",
            color="r",
            label="Random Precision",
        )
        plt.legend(loc="lower left")
        plt.title("Precision-Recall Curve")
        plt.xlabel("Recall")
        plt.ylabel("Precision")
        plt.xscale("log")
        plt.grid(True)
        plt.show()
    return metrics


def compute_epr(clf, X_test: np.ndarray, y_test: np.ndarray) -> float:
    """
    compute_epr computes the Expected Precision Recall (EPR) metric for the given classifier, test data, and true labels.

    Args:
        clf (sklearn.base.ClassifierMixin): The classifier to evaluate.
        X_test (numpy.ndarray): The test data features.
        y_test (numpy.ndarray): The true labels for the test data.

    Returns:
        float: The computed Expected Precision Recall (EPR) metric.
    """
    prb = clf.predict_proba(X_test)[:, 1]

    K = sum(y_test)
    # get only the top-K elems from prb
    pred = np.zeros(prb.shape)
    pred[np.argsort(prb)[-int(K) :]] = 1

    true_positive = np.sum(pred[y_test == 1] == 1)
    false_positive = np.sum(pred[y_test == 0] == 1)
    false_negative = np.sum(pred[y_test == 1] == 0)
    true_negative = np.sum(pred[y_test == 0] == 0)
    odds_ratio = (true_positive * true_negative) / (false_positive * false_negative)
    return odds_ratio


def unnormalize(df, is_root=False):
    # rescale (logp1 or sqroot)
    df = df**2 if is_root else (2**df) - 1
    r = np.array([i[i != 0].min() for k, i in df.iterrows()])
    df = (df.T / r).T.round()
    return df


def load_genes(organisms: Union[str, list] = "NCBITaxon:9606"):  # "NCBITaxon:10090",
    """
    load_genes loads the genes for the given organisms.

    Args:
        organisms (Union[str, list], optional): The organism(s) to load genes for. Can be a single organism string or a list of organism strings. Defaults to "NCBITaxon:9606".

    Returns:
        pd.DataFrame: A DataFrame containing gene information for the specified organisms, including columns for gene symbols, mitochondrial genes, ribosomal genes, hemoglobin genes, and organism.
    """
    try:
        import bionty as bt
    except ImportError:
        raise ImportError(
            "bionty is not installed. Please install it with pip install bionty \
            you will also need to populate its genes, have a look at jkobject/scdataloader's package \
            and its populate_ontology function"
        )
    organismdf = []
    if type(organisms) is str:
        organisms = [organisms]
    for organism in organisms:
        genesdf = bt.Gene.filter(
            organism_id=bt.Organism.filter(ontology_id=organism).first().id
        ).df()
        genesdf = genesdf[~genesdf["public_source_id"].isna()]
        genesdf = genesdf.drop_duplicates(subset="ensembl_gene_id")
        genesdf = genesdf.set_index("ensembl_gene_id").sort_index()
        # mitochondrial genes
        genesdf["mt"] = genesdf.symbol.astype(str).str.startswith("MT-")
        # ribosomal genes
        genesdf["ribo"] = genesdf.symbol.astype(str).str.startswith(("RPS", "RPL"))
        # hemoglobin genes.
        genesdf["hb"] = genesdf.symbol.astype(str).str.contains(("^HB[^(P)]"))
        genesdf["organism"] = organism
        organismdf.append(genesdf)
    organismdf = pd.concat(organismdf)
    organismdf.drop(
        columns=["source_id", "run_id", "created_by_id", "updated_at", "stable_id"],
        inplace=True,
    )
    return organismdf
