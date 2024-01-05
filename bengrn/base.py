"""
bengrn base module.

This is the principal module of the bengrn project.
here you put your main classes and objects.

Be creative! do whatever you want!

If you want to replace this with a Flask application run:

    $ make init

and then choose `flask` as template.
"""
from grnndata import GRNAnnData
from anndata import AnnData
from typing import Optional

# example constant variable
NAME = "bengrn"


class BenGRN:
    def __init__(self, adata: GRNAnnData, full_dataset: Optional[AnnData] = None):
        self.adata = adata
        self.full_dataset = full_dataset

    def do_tests(self, ground_truth: Optional[AnnData] = None):
        print("I'm doing tests!")


def scalefreeness():
    pass


def BeeLine():
    pass


def toGT():
    pass


def get_corr():
    pass


def get_genie3():
    pass


def load_GT_db(name="regulonDB", cell_type=None, organism="human"):
    pass
