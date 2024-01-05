# bengrn

[![codecov](https://codecov.io/gh/jkobject/benGRN/branch/main/graph/badge.svg?token=benGRN_token_here)](https://codecov.io/gh/jkobject/benGRN)
[![CI](https://github.com/jkobject/benGRN/actions/workflows/main.yml/badge.svg)](https://github.com/jkobject/benGRN/actions/workflows/main.yml)

Awesome Benchmark of Gene Regulatory Networks created by @jkobject

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

BaseClass(grndata).do_tests()
#or
some_test_function(grndata)
```

see more in the notebooks in the docs folder or in the [documentation](https://jkobject.github.io/benGRN/)

## Development

Read the [CONTRIBUTING.md](CONTRIBUTING.md) file.
