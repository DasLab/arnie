# Pseudoknots

[Pseudoknots](https://en.wikipedia.org/wiki/Pseudoknot) are a more complex form of secondary structure. 

<img src='/gh-docs/assets/pseudoknot.png' alt="example of a pseudoknot"></img>

Unpaired bases in a loop structure may pair with nucleotides elsewhere in the RNA sequence. This type of pairing is impossible to represent with the `(`, `)`, `.` and characters in traditional dot bracket notation, so we introduce new characters to represent various levels of pseudoknot pairings. In order, Arnie uses `[`, `{`, `<`, and lower case alphabet characters (`abc...`) to represent opening pairs, and `]`, `}`, `>`, and upper case alphabet characters (`ABC...`) to represent closing pairs. 

Here is an example pseudoknotted structure in dot bracket notation utilizing the expanded character set `...(((..[[[.(((...))))))]]]...`.

Many traditional structure prediction algorithms struggle with predicting pseudoknot structures, but there are a variety of approaches that can predict these complex folds. Arnie provides two main functions to predict pseudoknots: `pk_predict` and `pk_predict_from_bpp`.

## pk_predict
`pk_predict` takes an input RNA sequence string and returns a predicted secondary structure string in dot bracket notation that may include pseudoknots. It's very similar to the `mfe` function, but supports a different set of predictor packages that focus on pseudoknot prediction. 

**Args:**
```
  seq (str): nucleic acid sequence, required
  predictor (str): the folding library to use
  ipknot options:
    model: one of ["LinearPartition-C","LinearPartition-V","Boltzmann","ViennaRNA","CONTRAfold","NUPACK"]
    t1: probability threshold level 1 
    t2: probability threshold level 2
    refinement: number of times for refinment

  hotknots options:
    model: one of ["CC","RE","DP"]
    param: one of ["parameters_CC06.txt","parameters_CC09.txt","parameters_DP03.txt","parameters_DP09.txt"]

  spotrna options:
    cpu: number cpu threads
```

**Returns:**
```
  A string in dot bracket notation representing the predicted secondary structure of the provided sequence, potentially including pseudoknots.
```

**Example:** 
```
pk_predict("GUAUCAAAAAAGAUACGCCGUAUGCUAAUAUGUAUCUAUACUUGCUCUACAGGUUGAG", "knotty")

'..........(((((([[[[[[.[[...[[[))))))]]]...]]..]]].]]]....'
```

**Supported packages:**
- `hotknots`
- `ipknot`
- `knotty`
- `spotrna`
- `spotrna2`
- `e2efold`
- `pknots`
- `nupack`

## pk_predict_from_bpp
`pk_predict_from_bpp` takes a different approach to pseudoknot prediction. Rather than use dedicated pseudoknot prediction packages, `pk_predict_from_bpp` uses post-processing algorithms that can predict likely pseudoknots based on a sequence's predicted base pair probability matrix. This allows us to examine sequences for predicted pseudoknots with traditional predictive models that don't support pseudoknots by default. 

`pk_predict_from_bpp` provides two processing algorithms, [`threshknot`](https://arxiv.org/abs/1912.12796) and [`hungarian`](https://en.wikipedia.org/wiki/Hungarian_algorithm).

**Args:**
```
  bpp (array): base pair probability matrix, required
  heuristic (str): the pk prediction algorithm to use; either "hungarian" or "threshknot"
  threshknot options:
    theta
    max_iter
    allowed_buldge_len
    min_len_helix

  hungarian options:
    add_p_unpaired
    theta (aka prob_to_0_threshold_post)
    prob_to_0_threshold_prior
    prob_to_1_threshold_prior
    exp
    sigmoid_slope_factor
    ln
    allowed_buldge_len
    min_len_helix
```

**Returns:**
```
  A string in dot bracket notation representing the predicted secondary structure of the provided sequence, potentially including pseudoknots.
```

**Example:** 
```
bpps = bpps("GUAUCAAAAAAGAUACGCCGUAUGCUAAUAUGUAGGCGCUAUACUUGCUCUACACCGGCGGUUGAG", package="eternafold")
pk_predict_bpp(bpps)

'(((((......)))))..........................................'
```

**Supported packages:**
- `eternafold`
- `contrafold`
- `vienna`
- `nupack`
- `rnasoft`
- `rnastructure`
- `vfold`

