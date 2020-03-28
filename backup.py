def combine_constraint(ss_list):
    """ Combines a list of dot-bracket notation secondary structure 

    returns
        combined secondary structure constraint
        or False if not possible (or results in a pesudo knot)
    """
    num_ss = len(ss_list)
    if num_ss == 1:
        return ss_list[0]
    if len(np.unique([len(x) for x in ss_list])) > 1:
        raise ValueError('cannot combine ss constraints')
    dbn_string = ['.']*len(ss_list[0])
    merge_order = [] # order that the secondary structures are merged
    for idx, row in enumerate(zip(*ss_list)):
        if row.count('.') < num_ss-1:
            # Multiple secondary structure constraints cannot be merged
            return False
        if row.count('.') == num_ss:
            # If no constraints move to next position
            continue
        ss_idx, c = [[x,y] for x,y in enumerate(row) if y != '.'][0]
        dbn_string[idx] = c
        merge_order.append(ss_idx)
    merge_order = np.array(merge_order)
    
    # Check that aptamer secondary structures are not intermixed (NO PSUEDO KNOTS)
    # because packages cannot handle pseudo knots
    ss_idx_set = set(range(num_ss))
    cur_len = num_ss
    prev_len = num_ss
    while len(ss_idx_set):
        for ss_idx in ss_idx_set.copy(): # remove continuous stretches of numbers
            idx = np.where(merge_order == ss_idx)[0]
            dx = np.diff(idx)
            if len(idx) == 1 or np.all(dx == 1):
                merge_order = np.delete(merge_order, np.s_[idx[0]:idx[-1]+1])
                ss_idx_set.remove(ss_idx)
        cur_len = len(ss_idx_set)
        if cur_len == prev_len: # Unable to resolve. Must be psuedoknot
            return False
        prev_len = cur_len
    
    dbn_string = ''.join(dbn_string)
    if '()' in dbn_string: # () is impossible
        return False 
    else:
        return dbn_string

def prune_combo_list(combo_list, apt_idx_list, N_frag):
    """ Helper function for write_constraints 
    Prunes all possible combinations of the aptamer fragments
    """
    final_idx_list = []
    for idx, combo in enumerate(combo_list):
        temp = apt_idx_list[combo,:]
        start_idx_list = temp[:,1]
        ss_len_list = temp[:,3]
        if N_frag == 1: # If <= 2 fragments than always in order
            final_idx_list.append(idx)
        elif N_frag == 2:
            dx = np.diff(start_idx_list)
            if dx != 0:
                final_idx_list.append(idx)
        else:
            # first check if aptamer is in the correct order in the sequence
            dx = np.diff(start_idx_list)
            if np.all(dx < 0) or np.all(dx > 0):
                # second check if aptamer fragments don't overlap
                # each element in dx must be >= length of ss
                if dx[0] > 0: # fragments in increasing order
                    ss_len_list = ss_len_list[:-1]
                else: # fragments in decreasing order
                    ss_len_list = ss_len_list[1:]
                
                # Check if there is enough bases to fit the aptamer fragments
                if np.all(np.abs(dx) >= ss_len_list):
                    final_idx_list.append(idx)
    final_combo_list = combo_list[final_idx_list]
    return np.array(final_combo_list)

def combo_list_to_dbn_list(seq, final_combo_list, apt_idx_list, apt_ss_list):
    """ Helper function for write_constraints 
    Converts the combination of aptamer fragments to dbn strings
    """
    dbn_string_list = []
    for combo in final_combo_list:
        dbn_string=['.']*len(seq)
        temp = apt_idx_list[combo,:][:,1:3]
        
        # Fill in the dbn string with the aptamer ss
        for (start, finish), apt_ss in zip(temp, apt_ss_list):
            dbn_string[start:finish] = list(apt_ss)
        dbn_string = ''.join(dbn_string)
        
        # Now look at the 1st parenthesis to see if aptamer is backwards
        if dbn_string.find('(') > dbn_string.find(')'): # aptamer is reversed
            dbn_string = dbn_string.replace('(','?')
            dbn_string = dbn_string.replace(')','(')
            dbn_string = dbn_string.replace('?',')')
        if '()' not in dbn_string and len(dbn_string) > 0: # Check for impossible basepairs
            dbn_string_list.append(dbn_string)
    return dbn_string_list

def write_combo_constraints(seq, raw_apt_seq, raw_apt_ss, verbose=False):
    """ Given a sequence it will give all possible secondary constraints of the aptamer 

    Args:
        seq: RNA sequence
        raw_apt_seq: list of aptamer sequences
            e.g. ['CCC+GGG','AAUU', ...]
            + denotes splitable aptamer
        raw_apt_ss: list of aptamer secondary structure
            e.g. ['(((+)))','(xx)', ...]
            + denotes splitable aptamer
            x denotes unpaired base
            . denotes wildcard (can be anything)
        verbose: to be verbose
    Returns
        list of list of secondary structure constraints
            the 1st list corresponds to the 1st aptamer
            the 2nd list corresponds to the 2nd aptamer
            etc.
        Will contains ss constraints for multiple copies as well
        Will NOT return any constraints that contain ()
    """
    if verbose:
        print('Writing constraints')
        print(seq)
        print(raw_apt_seq)
        print(raw_apt_ss)
    # split everything by +
    apt_seq_list = raw_apt_seq.split('+')
    apt_ss_list = raw_apt_ss.split('+')
    if len(apt_seq_list) != len(apt_ss_list):
        raise ValueError('Missing + in aptamer sequence or secondary structure')

    # Iterate through each aptamer fragment and save its idx,
    apt_idx_list = []
    for idx, (apt_seq, apt_ss)  in enumerate(zip(apt_seq_list, apt_ss_list)):
        if seq.find(apt_seq) == -1:
            raise ValueError("Aptamer segment {} not found".format(idx+1))
        if len(apt_seq) != len(apt_ss):
            raise ValueError("Mismatch between aptamer sequence and aptamer secondary structure")

        # Save all locations of each fragment
        for m in re.finditer(apt_seq, seq):
            #Note: cannot get overlapping segments
            span = m.span()
            start = span[0]
            finish = span[1]
            ss_length = len(apt_seq)
            apt_idx_list.append([idx, start, finish, ss_length])
    apt_idx_list = np.array(apt_idx_list)

    # Now combine aptamer fragments into full secondary structure constraints
    N_frag = len(apt_ss_list) # Number of fragments to stitch together

    # Get a list of all possible combination (each combination is a list as well)
    temp = np.array([np.where(apt_idx_list[:,0] == idx)[0] for idx in range(N_frag)])
    if len(temp) > 1:
        combo = np.meshgrid(*temp)
        combo_list = np.array(combo).T.reshape(-1,N_frag) # reformat the combinations
    else:
        combo_list = np.array([[x] for x in temp[0]])

    # Check each combination to make sure it's feasible, if not remove
    final_combo_list = prune_combo_list(combo_list, apt_idx_list, N_frag)
    
    # Convert each combination into a dbn_string
    dbn_list = combo_list_to_dbn_list(seq, final_combo_list, apt_idx_list, apt_ss_list)
    
    # List will be empty if nothing can be in order (only for >= 3 fragments)
    return dbn_list