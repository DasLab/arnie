import numpy as np
from arnie.bpps import bpps

def threshknot_util(sequence, package='vienna_2', theta=0):
    '''
    Inputs:
    sequence: RNA sequence
    package: folding package to use
    
    Set theta = 0 to not filter base pairs as in ThreshKnot.
    
    Returns: N x N matrix of base pair probabilities. Nonzero entries represent base pairs
    predicted in final (possibly pseudoknotted) structure.
    Probabilities are their associated probability (obvs).
    '''
    
    bp_matrix = bpps(sequence, package=package)
    
    # if desired, filter base pair probabilities below a cutoff
    bp_matrix[np.where(bp_matrix <= theta)] = 0
    output = np.zeros([len(sequence),len(sequence)])
    
    # ProbKnot heuristic part 1: get all base pairs where p(ij) == p_max(i)
    output[np.where(bp_matrix == np.max(bp_matrix,axis=0))] = 1
    
    # ProbKnot heuristic part 2: get all base pairs where p(ij) == p_max(j)
    array_of_bps = np.clip(output+np.transpose(output)-1,0,1)
    
    # setting all bp probabilities not corresponding to a final selected base pair to zero
    bp_matrix[np.where(array_of_bps == 0)] = 0
    
    return bp_matrix
