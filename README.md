# arnie
Python API to compute RNA energetics and do structure prediction from all available packages.

Das Lab, 2019

Hannah Wayment-Steele, with inspiration from MWu's `wuami_tools`

## Organization:

`examples`: example jupyter notebooks with usage.

`scripts`: scripts for processing sequences in batch.

`test`: unit tests (still in work)

`mea`: code for computing Maximum Expected Accuracy structures.

## Usage:

Start by creating a file that points to your builds of all the structure prediction packages you wish to use.  An example file is provided in `example_arnie_file.txt`.  Create a variable in your .bashrc for this:

```
export ARNIEFILE="/path/to/arnie/<my_file.txt>"
```
NB: this file is technically yaml format, but isn't read in by yaml.

See `examples/start_here.ipynb` for example syntax. In brief, comparing across packages is simple. For computing base pairing probability matrices:

```
from arnie.bpps import bpps
%pylab inline
example_seq = 'GGGGAAAACCCC'

bpps={}

for pkg in ['vienna','nupack','RNAstructure','contrafold','RNAsoft']:
    bpps[package] = bpps(example_seq, package=pkg)
    
imshow(bpps['vienna'])
```

## Riboswitch fold change:

See `examples/riboswitch_fold_change.ipynb` for example k_d prediction and fold change prediction.

## Coming soon

Help for compiling packages
