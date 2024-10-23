import os

import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from grnndata import GRNAnnData
from scipy.sparse import csr_matrix

from bengrn.base import (
    NAME,
    BenGRN,
    compute_epr,
    compute_genie3,
    get_GT_db,
    get_perturb_gt,
    get_sroy_gt,
    train_classifier,
)


def test_base():
    assert NAME == "bengrn"
    adata = sc.read_h5ad(os.path.join(os.path.dirname(__file__), "test.h5ad"))
    random_matrix = np.random.rand(1000, 1000)
    adata = adata[:, np.argsort(-adata.X.sum(axis=0)).tolist()[0][:1000]]
    random_mask = np.random.choice([0, 1], size=random_matrix.shape, p=[0.8, 0.2])
    sparse_random_matrix = csr_matrix(random_matrix * random_mask)
    grn = GRNAnnData(adata.copy(), grn=sparse_random_matrix)
    grn.var.index = grn.var.symbol.astype(str)
    _ = BenGRN(grn, doplot=False).scprint_benchmark()

    # Test get_sroy_gt function
    sroy_gt = get_sroy_gt(get="liu")
    assert isinstance(
        sroy_gt, GRNAnnData
    ), "get_sroy_gt should return a GRNAnnData object"

    # Test get_perturb_gt function
    perturb_gt = get_perturb_gt()
    assert isinstance(
        perturb_gt, GRNAnnData
    ), "get_perturb_gt should return a GRNAnnData object"

    # Test compute_genie3 function
    genie3_result = compute_genie3(adata[:, :100], ntrees=10, nthreads=1)
    assert isinstance(
        genie3_result, GRNAnnData
    ), "compute_genie3 should return a GRNAnnData object"

    # Test train_classifier function
    random_matrix = np.random.rand(4, 10000).reshape(100, 100, 4)
    subgrn = grn[:, :100]
    subgrn.varp["GRN"] = random_matrix
    classifier, metrics, clf = train_classifier(subgrn)
