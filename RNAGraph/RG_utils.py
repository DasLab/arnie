from collections import Counter
import numpy as np
from copy import copy

#Copied from pyFold-1D (H Wayment-Steele)

def find_all(s, ch):
    if isinstance(ch, list):
        tmp_list=[]
        for c in ch:
            tmp_list.extend([i for i, ltr in enumerate(s) if ltr == c])
        return tmp_list
    elif isinstance(ch, str):
        return [i for i, ltr in enumerate(s) if ltr == ch]

def convert_structure_to_bps(secstruct):

    bps = []

    #Find other delimiters
    other_delimiters = [k for k in Counter(secstruct).keys() if k not in ['.','(',')','[',']','{','}']]

    for delim in other_delimiters:
        pos = find_all(secstruct, delim)
        assert(len(pos) % 2 == 0)

        N = int(len(pos)/2)
        i,j = pos[:N], pos[N:]

        if N > 1:
            assert(all(np.diff(i)==1))
            assert(all(np.diff(j)==1))

        for ind in range(N):
            bps.append([i[ind],j[-1-ind]])

    left_delimiters = ['(','{','[']
    right_delimiters = [')','}',']']

    for (left_delim, right_delim) in list(zip(left_delimiters, right_delimiters)):

        left_list = []
        for i, char in enumerate(secstruct):
            if char == left_delim:
                left_list.append(i)

            elif char == right_delim:
                bps.append([left_list[-1],i])
                left_list = left_list[:-1]

        assert len(left_list)==0

    return bps

def parse_stems_from_bps(bps, debug=False):
    #Creates list of stems, where each stem is a list of base pairs.
    # i.e. '((.))' -> [[0,4],[1,3]]
    # '((.)).((.))' -> [[[0,4],[1,3]],[[6,10],[7,9]]]
    
    if debug: print(bps)

    if len(bps) == 0:
        stems = []
    else:
        nres = np.max(bps)
        stems = []

        while len(bps) > 0:
            bp = bps[0]
            bps = bps[1:]

            stem = [bp]

            if debug: print('stem init', stem)
            bp_next = copy(bp)
            if debug: print('bp_next', bp_next)

            # Check outward
            for i in list(reversed(range(bp[0]))):
                bp_next = [copy(bp_next)[0]-1,copy(bp_next)[1]+1]

                if debug: print('next_out', bp_next)
                if debug: print('bps here', bps)
                if len(bps) > 0:
                    gp = find_all([x[0] for x in bps], [bp_next[0]])
                    if debug: print('gp', gp)
                    if len(gp)>0:
                        if bps[gp[0]][1] == bp_next[1]: # found an extension
                            if debug: print('r')
                            stem.append(copy(bp_next))
                            del bps[gp[0]] # take out of bp list
                    else:
                        break

            bp_next = copy(bp)

            #Check inward
            for i in range(bp[0],nres+1):

                bp_next[0] = copy(bp_next)[0]+1
                bp_next[1] = copy(bp_next)[1]-1

                if debug: print('next_in', bp_next)
                if len(bps) > 0:
                    gp = find_all([x[0] for x in bps], [bp_next[0]])
                    if len(gp)>0:
                        if bps[gp[0]][1] == bp_next[1]: # found an extension
                            if debug: print('h')
                            stem = [copy(bp_next)]+copy(stem)
                            del bps[gp[0]] # take out of bp list
                        else:
                            break
            stems.append(stem)
            if debug: print('stem', stem)
    if debug: print('stems', stems)
    return stems

def parse_out_chainbreak(secstruct):
    secstruct_new = []
    is_chainbreak = []

    for char in secstruct:
        if char in [',','+',' ','&']:
            if len(is_chainbreak)>0:
                is_chainbreak[-1] = 1
        else:
            secstruct_new.append(char)
            is_chainbreak.append(0)
    return is_chainbreak, ''.join(secstruct_new)


def get_stem_assignment(secstruct):
    '''Returns vector length N_beads, 0 if not in a stem, otherwise assigned stem 1 to max(N_stems)
    Note basically not switched to zero-indexing for python version, unlike partner syntax'''

    is_chainbreak, secstruct = parse_out_chainbreak(secstruct)

    if not secstruct[0].isdigit():
        bps = convert_structure_to_bps(secstruct)

    else: # assume partner vector was given
        partner = secstruct
        bps = []
        for i in range(len(partner)):
            if partner[i] > i:
                bps.append([i, partner[i]])

    stems = parse_stems_from_bps(bps)

    stem_assignment = np.zeros([len(secstruct)])

    for i, stem in enumerate(stems):
        for bp in stem:
            stem_assignment[bp[0]] = i+1
            stem_assignment[bp[1]] = i+1

    return stem_assignment

def get_pairmap(secstruct):
    """
    (lifted from draw_rna)
    generates a list containing pair mappings

    args:
    secstruct contains secondary structure string

    returns:
    list with pair mappings
    """
    pair_stack = []
    end_stack = []
    pairs_array = []
    i_range = list(range(0,len(secstruct)))

    # initialize all values to -1, meaning no pair
    for ii in i_range:
        pairs_array.append(-1)

    # assign pairs based on secstruct
    for ii in i_range:
        if(secstruct[ii] == "("):
            pair_stack.append(ii)
        elif(secstruct[ii] == ")"):
            if not pair_stack:
                end_stack.append(ii)
            else:
                index = pair_stack.pop()
                pairs_array[index] = ii
                pairs_array[ii] = index
    if len(pair_stack) == len(end_stack):
        n = len(pair_stack)
        for ii in range(n):
            pairs_array[pair_stack[ii]] = end_stack[-ii]
            pairs_array[end_stack[-ii]] = pair_stack[ii]
    else:
         print("ERROR: pairing incorrect %s" % secstruct)

    return pairs_array