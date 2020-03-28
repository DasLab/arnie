import numpy as np
import argparse, sys

def convert_dotbracket_to_matrix(s):
    m = np.zeros([len(s),len(s)])
    for char_set in [['(',')'], ['[',']'],['{','}']]:
        bp1=[]
        bp2=[]
        for i, char in enumerate(s):
            if char==char_set[0]:
                bp1.append(i)
            if char==char_set[1]:
                bp2.append(i)
        for i in list(reversed(bp1)):
            for j in bp2:
                if j > i:
                    m[i,j]=1.0
                    bp2.remove(j)
                    break
    return m

def convert_matrix_to_dotbracket(m):
    print('not implemented yet!')
    pass

def load_matrix_or_dbn(s):
    num_lines = sum(1 for line in open(s))

    if num_lines > 2: #heuristic here
        struct = np.loadtxt(s) # load as base pair matrix
        assert struct.shape[0] == struct.shape[1]
    else:
        try: # load as dot-bracket string

            dbn_struct = open(s,'r').read().rstrip()
            
            struct = convert_dotbracket_to_matrix(dbn_struct)
        except:
            raise ValueError('Unable to parse structure %s' % s)
    return struct

def score_ground_truth(pred_matrix, true_matrix):
    '''Score a predicted structure against a true structure,
     input as NxN base pair matrix (takes top triangle).'''

    N = pred_matrix.shape[0]
    #print('pred',pred_matrix.shape, 'true', true_matrix.shape)
    assert pred_matrix.shape[1] == N
    assert true_matrix.shape[0] == N
    assert true_matrix.shape[1] == N

    true = true_matrix[np.triu_indices(N)]
    pred = pred_matrix[np.triu_indices(N)]

    TP, FP, cFP, TN, FN = 0, 0, 0, 0, 0

    for i in range(len(true)):
        if true[i] == 1:
            if pred[i] == 1: 
                TP += 1
            else:
                FN += 1
        elif true[i] == 0:
            if pred[i] == 0: 
                TN += 1
            else: 
                FP += 1
                #check for compatible false positive
                a,b = np.triu_indices(N)
                if np.sum(true_matrix,axis=0)[a[i]]+ np.sum(true_matrix,axis=0)[b[i]]==0:
                   cFP +=1

    # cFP = 0 #for debugging

    #print('TP', TP, 'TN', TN, 'FP', FP, 'FN', FN, 'cFP', cFP)

    if TP + FN == 0:
        sen = 1
    else:
        sen = TP/(TP + FN)

    if TP + FP - cFP == 0:
        ppv = 1
    else:
        ppv = TP/(TP + FP - cFP)

    mcc_num = (TP*TN - (FP - cFP)*FN)
    mcc_denom = np.sqrt((TP + FP - cFP)*(TP + FN)*(TN + FP - cFP)*(TN + FN))

    if  mcc_denom == 0:
        mcc = mcc_num
    else:
        mcc = mcc_num/mcc_denom

    if ppv + sen == 0:
        fscore = 0
    else:
        fscore = 2*ppv*sen/(ppv+sen)

    return sen, ppv, mcc, fscore, N
