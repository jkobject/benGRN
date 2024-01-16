# bengrn

[![codecov](https://codecov.io/gh/jkobject/benGRN/branch/main/graph/badge.svg?token=benGRN_token_here)](https://codecov.io/gh/jkobject/benGRN)
[![CI](https://github.com/jkobject/benGRN/actions/workflows/main.yml/badge.svg)](https://github.com/jkobject/benGRN/actions/workflows/main.yml)

Awesome Benchmark of Gene Regulatory Networks created by @jkobject

The package is supposed to work with [GRnnData](https://cantinilab.github.io/GRnnData/)
It can run Genie3 & pyscenic on your data as a comparison

It has 4 different types of key ground truth data to compare your GRN to:
- sushmita roy's ChIP+Perturb ground truth
- collectri's literature curated ground truth
- dorothea's literature curated ground truth
- tf2gene's chip curated ground truth

You can find the documentation [here](https://jkobject.com/benGRN/)

## Install it from PyPI

```bash
pip install bengrn
```

## Usage

```py
from bengrn import BenGRN
from bengrn import some_test_function

# a GRN in grnndata formart
grndata

BenGRN(grndata).do_tests()
#or
some_test_function(grndata)
```

see more in the notebooks in the docs folder or in the [documentation](https://jkobject.com/benGRN/)

## Development

Read the [CONTRIBUTING.md](CONTRIBUTING.md) file.
