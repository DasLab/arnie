# arnie
Arnie is a Python API to compute RNA energetics and do structure prediction across multiple secondary structure packages.

## Install
`arnie` is [available on PyPI](https://pypi.org/project/arnie/).

`pip install arnie`

## Simple Setup
Arnie works by delegating calls for structure predictions to various RNA prediction libraries. To use arnie we need to have these libraries installed, and we need to point to these their installed locations with environment variables. Here we will use [Eternafold](https://github.com/eternagame/Eternafold) which is simple to install via [Bioconda](https://bioconda.github.io/recipes/eternafold/README.html). This example assumes you have conda installed already; see the full [setup page](/setup/environment.md) for more details about setting up an arnie environment.


```
conda install -c bioconda eternafold
export eternafold_PATH=/path/to/installed/location
```

## Usage:

See the [usage docs](/usage/structure_prediction) for example syntax. In brief, comparing across packages is simple. For computing base pairing probability matrices:

```
from arnie.bpps import bpps

bpps_dict = {}
my_sequence = 'CGCUGUCUGUACUUGUAUCAGUACACUGACGAGUCCCUAAAGGACGAAACAGCG'

for pkg in ['vienna','nupack','RNAstructure','contrafold','RNAsoft']:
    bpps_dict[pkg] = bpps(my_sequence, package=pkg)
```

(c) 2024 [Das Lab](https://daslab.stanford.edu/), Leland Stanford Jr University