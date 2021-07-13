import numpy as np
import argparse, sys

def convert_dotbracket_to_matrix(s):
    m = np.zeros([len(s),len(s)])
    for char_set in [['(',')'], ['[',']'],['{','}'],['<','>']]:
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
    bp_list = convert_matrix_to_bp_list(m)
    return convert_bp_list_to_dotbracket(bp_list,len(m))

def convert_matrix_to_bp_list(m):
    bp_list = [] # convert adjacency matrix to adjacency list
    for i,row in enumerate(m):
        for j,is_bp in enumerate(row[i+1:]):
            if is_bp:
                bp_list.append([i,i+1+j])
    return bp_list


def convert_bp_list_to_dotbracket(bp_list,seq_len):
    dotbracket = "."*seq_len
    # group into bps that are not intertwined and can use same brackets!
    groups = group_into_non_conflicting_bp_(bp_list)

    # all bp that are not intertwined get (), but all others are
    # groups to be nonconflicting and then asigned (), [], {}, <> by group
    chars_set = [("(",")"),("(",")"),("[","]"),("{","}"),("<",">")]
    if len(groups) > len(chars_set):
        print("WARNING: PK too complex, not enough brackets to represent it.")

    for group,chars in zip(groups,chars_set):
        for bp in group:
            dotbracket = dotbracket[:bp[0]] + chars[0] + dotbracket[bp[0]+1:bp[1]] + chars[1] + dotbracket[bp[1]+1:]
    return dotbracket


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


def group_into_non_conflicting_bp_(bp_list):
    ''' given a conflict list from get_list_bp_conflicts_, group basepairs into groups that do not conflict

    Args
        conflict_list: list of pairs of base_pairs that are intertwined basepairs

    Returns:
        groups of baspairs that are not intertwined
    '''
    conflict_list = get_list_bp_conflicts_(bp_list)

    non_redudant_bp_list = get_non_redudant_bp_list_(conflict_list)
    bp_with_no_conflict = [bp for bp in bp_list if bp not in non_redudant_bp_list]
    groups = [bp_with_no_conflict]
    while non_redudant_bp_list != []:
        current_bp = non_redudant_bp_list[0]
        current_bp_conflicts = []
        for conflict in conflict_list:
            if current_bp == conflict[0]:
                current_bp_conflicts.append(conflict[1])
            elif current_bp == conflict[1]:
                current_bp_conflicts.append(conflict[0])
        group = [bp for bp in non_redudant_bp_list if bp not in current_bp_conflicts]
        groups.append(group)
        non_redudant_bp_list = current_bp_conflicts
        conflict_list = [conflict for conflict in conflict_list if conflict[0] not in group and conflict[1] not in group]
    return groups


def get_list_bp_conflicts_(bp_list):
    '''given a bp_list gives the list of conflicts bp-s which indicate PK structure
    Args:
        bp_list: of list of base pairs where the base pairs are list of indeces of the bp in increasing order (bp[0]<bp[1])
    returns:
        List of conflicting basepairs, where conflicting is pairs of base pairs that are intertwined.
    '''
    if len(bp_list) <= 1:
        return []
    else:
        current_bp = bp_list[0]
        conflicts = []
        for bp in bp_list[1:]:
            if (bp[0] < current_bp[1] and current_bp[1] < bp[1]):
                conflicts.append([current_bp,bp])
        return conflicts + get_list_bp_conflicts_(bp_list[1:])

def get_non_redudant_bp_list_(conflict_list):
    ''' given a conflict list get the list of nonredundant basepairs this list has

    Args:
        conflict_list: list of pairs of base_pairs that are intertwined basepairs
    returns:
        list of basepairs in conflict list without repeats
    '''
    non_redudant_bp_list = []
    for conflict in conflict_list:
        if conflict[0] not in non_redudant_bp_list:
            non_redudant_bp_list.append(conflict[0])
        if conflict[1] not in non_redudant_bp_list:
            non_redudant_bp_list.append(conflict[1])
    return non_redudant_bp_list
