# arnie
Python API to compute RNA energetics and do structure prediction across multiple secondary structure packages.

Currently supported:

Vienna^1 [https://www.tbi.univie.ac.at/RNA/#download]

NUPACK^2 [http://www.nupack.org/downloads]

RNAstructure^3 [https://rna.urmc.rochester.edu/RNAstructure.html]

RNAsoft^4 [http://www.rnasoft.ca/download.html]

CONTRAfold^5 [http://contra.stanford.edu/contrafold/]

(c) 2020 Leland Stanford Jr University

Authors:
Hannah Wayment-Steele

## Organization:

`notebooks`: example jupyter notebooks with usage.

`scripts`: scripts for processing sequences in batch.

`parameter_files`: dir of various parameter files for packages, put here out of convenience.

`test`: unit tests (still in work)

`mea`: code for computing Maximum Expected Accuracy structures.

`RNAGraph`: Code to process secondary structures as graph objects (here until a better repo location becomes clear)

## Usage:

Start by creating a file that points to your builds of all the structure prediction packages you wish to use.  An example file is provided in `example_arnie_file.txt`.  Create a variable in your .bashrc for this:

```
export ARNIEFILE="/path/to/arnie/<my_file.txt>"
```
NB: this file is technically yaml format, but isn't read in by yaml.

See `examples/start_here.ipynb` for example syntax. In brief, comparing across packages is simple. For computing base pairing probability matrices:

```
from arnie.bpps import bpps

bpps={}

for pkg in ['vienna','nupack','RNAstructure','contrafold','RNAsoft']:
    bpps[package] = bpps(my_sequence, package=pkg)
```
![](assets/example_base_pair_matrices.png)

Can also analyze as average base pairing per nucleotide:

![](assets/example_avg_bp_per_nucleotide.png)


## Coming soon

Pointers for compiling packages

## References

1. Lorenz, R. et al. ViennaRNA Package 2.0. Algorithms Mol Biol 6, 26 (2011).
2. Zadeh, J.N. et al. NUPACK: Analysis and design of nucleic acid systems. J Comput Chem 32, 170-173 (2011).
3. Reuter, J.S. & Mathews, D.H. RNAstructure: software for RNA secondary structure prediction and analysis. BMC Bioinformatics 11, 129 (2010).
4. Andronescu, M., Condon, A., Hoos, H.H., Mathews, D.H. & Murphy, K.P. in RNA, Vol. 16 2304-2318 (2010).
5. Do, C.B., Woods, D.A. & Batzoglou, S. CONTRAfold: RNA secondary structure prediction without physics-based models. Bioinformatics 22, e90-98 (2006).
