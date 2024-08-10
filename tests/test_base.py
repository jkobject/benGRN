import scanpy as sc
import numpy as np
import os
from scipy.sparse import csr_matrix

from grnndata import GRNAnnData
from bengrn.base import NAME, BenGRN


def test_base():
    assert NAME == "bengrn"
    adata = sc.read_h5ad(os.path.join(os.path.dirname(__file__), "test.h5ad"))
    random_matrix = np.random.rand(1000, 1000)
    adata = adata[:, np.argsort(-adata.X.sum(axis=0)).tolist()[0][:1000]]
    random_mask = np.random.choice([0, 1], size=random_matrix.shape, p=[0.8, 0.2])
    sparse_random_matrix = csr_matrix(random_matrix * random_mask)
    grn = GRNAnnData(adata.copy(), grn=sparse_random_matrix)
    grn.var.index = grn.var.symbol.astype(str)
    metrics = BenGRN(grn, doplot=False).scprint_benchmark()
    assert True, "BenGRN test passed"
