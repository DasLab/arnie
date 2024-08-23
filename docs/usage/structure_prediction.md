
## Structure Prediction

## MFE
The `mfe` function generates a "minimum free energy" structure prediction with the selected package. The minimum free energy prediction is the secondary structure calculated to have the lowest free energy value. In theory, the lower the free energy, the more likely the structure is to form. Not all predictors support free energy-based estimates (although many do). 

Note: `mfe` operates differently than [`mea`](#mea). That said, contrafold's default structure prediction is an MEA structure, not MFE. When using `mfe`, calling contrafold returns the default MEA structure unless the `--viterbi` flag is used, which will use the viterbi (MFE) algorithm in contrafold. 


**Args:**
```
  seq (str): nucleic acid sequence, required
  package (str): the folding library to use
  T (float): temperature (Celsius)
  constraint (str): structure constraints
  motif (str): argument to vienna motif 
  linear (bool): call LinearFold to estimate MFE in Vienna or Contrafold
  return_dG_MFE (bool): also return dG(MFE) (specific to linearfold)
  dangles (bool): dangles or not (specific to linearfold)
  noncanonical(bool): include noncanonical pairs or not (specific to contrafold, RNAstructure (Cyclefold))
  param_file(str): path to specific thermodynamic parameter file (specific to contrafold, eternafold)
  coaxial (bool): coaxial stacking or not (specific to rnastructure)
  viterbi (bool): use the viterbi algorithm for mfe calculation (specific to contrafold)
  pseudo (bool): if True, will predict pseudoknots
  shape_signal (list): list of normalized SHAPE reactivities, with negative values indicating no signal (specific to rnastructure)
  dms_signal (list): list of normalized DMS reactivities, with negative values indicating no signal (specific to rnastructure)
  shape_file (str): path to file containing shape_signal (specific to rnastructure)
  dms_file (str): path to file containing dms_signal (specific to rnastructure)
```

**Returns:**
```
  A string in dot bracket notation representing the calculated MFE structure of the provided sequence.
```

**Example:** 
```
mfe("GUAUCAAAAAAGAUAC")
'(((((......)))))'
```

**Supported packages:**
- `eternafold`
- `contrafold`
- `vienna`
- `rnastructure`
- `linearfold`

## BPPS
The `bpps` function calculates the "base pairing probability matrix" with the selected package. The base pairing probaility matrix is an NxN matrix (where N is the length of the RNA sequence), with the value of the `i,j` position representing the probability of the `i` nucleotide pairing with the `j` nucleotide.

**Args:**
```
  sequence (str): nucleic acid sequence, required
  package (str): the folding library to use
  constraint (str): structure constraint [vienna, contrafold, rnastructure]
  linear (bool): call LinearPartition to estimate Z in Vienna or Contrafold

  motif (str): argument to vienna motif
  pseudo (bool): (NUPACK only) include pseudoknot calculation
  dangles (bool): dangles or not, specifiable for vienna, nupack
  dna (bool): (NUPACK only) use SantaLucia 1998 parameters for DNA
  coaxial (bool): coaxial stacking or not, specifiable for rnastructure, vfold
  noncanonical(bool): include noncanonical pairs or not (for contrafold, RNAstructure (Cyclefold))
  beam size (int): Beam size for LinearPartition base pair calculation.
  DEBUG (bool): Output command-line calls to packages.
  threshknot (bool): calls threshknot to predict pseudoknots (for contrafold with LinearPartition)
  shape_signal (list): list of normalized SHAPE reactivities, with negative values indicating no signal (specific to rnastructure)
  dms_signal (list): list of normalized DMS reactivities, with negative values indicating no signal (specific to rnastructure)
  shape_file (str): path to file containing shape_signal (specific to rnastructure)
  dms_file (str): path to file containing dms_signal (specific to rnastructure)
```

**Returns:**
```
  array: NxN matrix of base pair probabilities
```

**Example:** 
```
bpps("GUAUCAAAAAAGAUAC")
array([[0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 3.77178e-04,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 4.39771e-04, 0.00000e+00,
        8.24776e-01],
       [0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 1.69534e-04, 2.01963e-04, 1.93469e-04, 2.05658e-04,
        2.01099e-04, 1.37709e-04, 5.21924e-04, 0.00000e+00, 8.42528e-01,
        0.00000e+00],
...
```

**Supported packages:**
- `eternafold`
- `contrafold`
- `vienna`
- `nupack`
- `rnasoft`
- `rnastructure`
- `vfold`