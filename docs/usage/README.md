# Using Arnie
Arnie's primary purpose is to simplify the process of making structure predictions for an RNA sequence with a variety of structure prediction libraries. 

## RNA structure
RNA molecules form complex three-dimensional shapes in nature. We represent these forms in three structure levels of increasing complexity. 

1. **Primary Structure**
The primary structure of an RNA molecule is the base identity of the various nucleotides that make up the molecule. This sequence string is typically written in the 5' to 3' direction.
Example: "AGUAUCAAAAAAGAUAC"

2. **Secondary Structure**
The secondary structure of an RNA molecule is the set of base paring interactions between nucleotides in an RNA molecule. There are multiple ways to computationally represent secondary structure, although arnie primarily uses two: the base pairing matrix and the dot bracket string.

  A ***base pairing matrix*** is an NxN matrix (where N is the length of the RNA sequence), with the value of the `i,j` position representing the probability of the `i` nucleotide pairing with the `j` nucleotide.

  A ***dot bracket string*** is a representation of secondary structure where `(` and `)` characters represent base pairs and `.` characters represent unpaired bases. For example, `((....))` in dot bracket notation indicates that the 1st nucleotide is paired with the 8th nucleotide, the 2nd nucleotide is paired with the 7th nucleotide, and the others are unpaired. More complex secondary structures can also be represented in dot bracket notation (see [Pseudoknots](usage/pseudoknots.md) for more details).

  Arnie provides [several methods to predict secondary structures](usage/structure_prediction.md).

3. **Tertiary Structure**
The tertiary structure is the three-dimensional structure of the RNA molecule, with each atom located in a 3D coordinate space. Arnie doesn't work with this level of structure.

## Examples
The easiest way to get started with arnie is trying out our example notebooks to explore the functionality arnie provides.

- [Basic Introduction / Install](https://github.com/daslab/arnie/blob/master/notebooks/IntroToArnie.ipynb)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/daslab/arnie/blob/master/notebooks/IntroToArnie.ipynb)
