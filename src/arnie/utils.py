import os
import re
import subprocess as sp
import random
import string
import numpy as np
import arnie


def complement_to_(string):
    base_pairing_dct = {'a': 'u', 'u': 'a', 'g': 'c', 'c': 'g', 't': 'a'}
    return ''.join(base_pairing_dct[x.lower()] for x in string[::-1])


def run_RNAPVmin(probing_signal, seq, LOC, DEBUG, tauSigmaRatio=1, shapeConversion='S'):
    reac_file = write_reactivity_file_vienna(probing_signal, seq)
    fname = write([seq])
    RNApvmin_command = ['%s/RNApvmin' % LOC, reac_file, '--shapeConversion=%s' %
                        shapeConversion, '--tauSigmaRatio=%f' % tauSigmaRatio]

    with open(fname) as f:
        if DEBUG:
            print(fname)
        if DEBUG:
            print(' '.join(RNApvmin_command))
        p = sp.Popen(RNApvmin_command, stdin=f, stdout=sp.PIPE, stderr=sp.PIPE)
    rnapvmin_stdout, rnapvmin_stderr = p.communicate()

    shape_file = filename()

    with open(shape_file, 'wb') as f:
        f.write(rnapvmin_stdout)

    if DEBUG:
        print('stdout')
        print(rnapvmin_stdout)
        print('stderr')
        print(rnapvmin_stderr)

    if p.returncode:
        raise Exception('RNApvmin failed: on %s\n%s' % (seq, stderr))

    return shape_file


###############################################################################
# File writing
###############################################################################


def convert_dbn_to_RNAstructure_input(seq, constraints, filename):
    assert(len(seq) == len(constraints))

    bp_list = convert_dotbracket_to_bp_dict(constraints)

    SS_list, pairs_list = [], []

    for i, (s, c) in enumerate(list(zip(seq, constraints))):
        if c == 'x':
            SS_list.append(i+1)
        elif c == '.':
            pass
        elif c == '(':
            pairs_list.append([i+1, bp_list[i]+1])
        elif c == ')':
            pass
        else:
            print('Error reading constraint string', c)

    with open('%s' % filename, 'w') as out:
        out.write('DS:\n-1\n')
        out.write('SS:\n')
        for x in SS_list:
            out.write('%d\n' % x)
        out.write('-1\n')
        out.write('Mod:\n-1\n')
        out.write('Pairs:\n')
        for x, y in pairs_list:
            out.write('%d %d\n' % (x, y))
        out.write('-1 -1')
        out.write('FMN:\n-1\nForbids:\n-1\n')


def write_vector_to_file(vector, outfile):
    for x in vector:
        outfile.write('%.3f\n' % x)
    return


def write_matrix_to_file(matrix, outfile):
    for x in matrix:
        outfile.write('\t'.join(['%.3f' % y for y in x])+'\n')
    return

def write_reactivity_file_RNAstructure(reactivities, fname=None):
    """ writes reactivities (either SHAPE or DMS) to file format used by RNAstructure

      ex:
      1 0.120768
      2 0.190510
      3 0.155776

    Args:
      reactivities (list): a list of normalized reactivity float data. 
      Negative numbers can be used to indicate no signal.
    """

    if fname is None:
        fname = '%s.SHAPE' % filename()
    with open(fname, 'w') as f:
        i = 1
        for reactivity in reactivities:
            if reactivity >= 0:
                f.write('%d %f\n' % (i, reactivity))
            i += 1
    return fname


def write_reactivity_file_contrafold(reactivities, sequence, fname=None):
    '''write reactivity (either SHAPE or DMS) to file format used by CONTRAfold.

      ex:
      1 U e1 0.120768
      2 G e1 0.190510
      3 U e1 0.155776

    Args:
      reactivities (list): a list of normalized reactivity float data. 
      sequence: RNA sequence
      Negative numbers can be used to indicate no signal.
    '''

    assert len(reactivities) == len(sequence)

    if fname is None:
        fname = '%s.bpseq' % filename()

    with open(fname, 'w') as f:
        i = 1
        for char, reactivity in list(zip(sequence, reactivities)):
            if reactivity > 0:
                f.write('%d %s e1 %.6f\n' % (i, char, reactivity))
            elif reactivity == 0:
                f.write('%d %s e1 0.0001\n' % (i, char))
            else:
                f.write('%d %s e1 0.0\n' % (i, char))

            i += 1
    return fname


def local_rand_filename(n=6):
    """generate random filename

    Args:
      n (int): number of characters
    """
    rand = ''.join([random.choice(string.ascii_lowercase) for _ in range(n)])
    return rand


def get_random_folder(n=6):
    """ generate randome foldername
    that does not exist in current folder"""
    out_folder = ''.join([random.choice(string.ascii_lowercase) for _ in range(n)])
    while os.path.isdir(out_folder):
        out_folder = ''.join([random.choice(string.ascii_lowercase) for _ in range(n)])
    tmpdir = load_package_locations()['TMP']
    out_folder = f'{tmpdir}/{out_folder}'
    return out_folder


def filename(n=6):
    """generate random filename

    Args:
      n (int): number of characters
    """
    rand = ''.join([random.choice(string.ascii_lowercase) for _ in range(n)])
    tmpdir = load_package_locations()['TMP']
    return '%s/%s' % (tmpdir, rand)


def write(lines, fname=None):
    """write lines to file

    Args:
      lines (list): line(s) to write to file
      fname (str): filename to write to
    """
    if fname is None:
        fname = '%s.in' % filename()
    with open(fname, 'w') as f:
        for line in lines:
            f.write('%s\n' % line)
    return fname

def convert_dbn_to_contrafold_input(seq, constraints, filename):
    constraint_list = write_constraint_string(seq, constraints)
    with open('%s' % filename, 'w') as out:
        for i in range(len(seq)):
            out.write('%d\t%s\t%d\n' % (i+1, seq[i], constraint_list[i]))


###############################################################################
# File reading
###############################################################################

def bpseq_to_bp_list(bpseq_file, header_length=1):
    ''' read a bpseq file into a bp_list
    assumes first line header_length is title 
    and starts reading only after that
    '''
    bps = open(bpseq_file).readlines()
    bps = bps[header_length:]
    bp_list = []
    for bp in bps:
        bp = bp.split()[::2]
        bp = [int(nt) for nt in bp]
        if bp[1] != 0 and bp[0] < bp[1]:
            bp_list.append([bp[0] - 1, bp[1] - 1])
    return bp_list


def ct_to_bp_list(ct_file, header_length=1):
    ''' read a ct file into a bp_list
    assumes first line header_length is title 
    and starts reading only after that
    '''
    bps = open(ct_file).readlines()
    bps = bps[header_length:]
    bp_list = []
    for bp in bps:
        bp = bp.split()[::4]
        if len(bp) != 0:
            bp = [int(nt) for nt in bp]
            if bp[1] != 0 and bp[0] < bp[1]:
                bp_list.append([bp[0] - 1, bp[1] - 1])
    return bp_list


def prob_to_bpp(prob_file):
    ''' read a .prob file and return a bpp
    '''
    return np.loadtxt(prob_file)


###############################################################################
# Package handling
###############################################################################

supported_packages = [
    "contrafold",
    "eternafold",
    "nupack",
    "RNAstructure",
    "RNAsoft",
    "vienna"
    # PK Predictors
    "e2efold",
    "hotknots", 
    "ipknot", 
    "knotty", 
    "pknots",
    "spotrna",
    "spotrna2",
]

def print_path_files():
    package_dct = load_package_locations()
    for key, v in package_dct.items():
        print(key, v)


def package_list():
    pkg_list = []
    package_dct = load_package_locations()
    for key, v in package_dct.items():
        if key != "TMP" and key.lower() != 'bprna':
            if not key.startswith('linear'):
                if key == 'eternafoldparams' and 'eternafold' not in pkg_list:
                    pkg_list.append('eternafold')
                else:
                    if v != "None":
                        pkg_list.append(key)
    return pkg_list


def load_package_locations(DEBUG=False):
    '''Set up  paths to RNA folding packages. Checks environment variables or a user-supplied file. 
    If using the file, specify this in your ~/.bashrc as $ARNIEFILE'''
    return_dct = {}
    package_path = os.path.dirname(arnie.__file__)

    if DEBUG:
        print(supported_packages)

    # Read from environment variables
    for package in supported_packages:
        envVar = f"{package.upper()}_PATH"
        # Nupack installation sets its own environment variables
        if package == "nupack":
            envVar = f"{package.upper()}HOME"
        path = os.environ.get(envVar)
        if path:
            return_dct[package] = path
            if DEBUG:
                print(f'{package}: {path}')

    # Read from arnie file as last resort
    if os.environ.get("ARNIEFILE"):
        if DEBUG:
            print('Reading Arnie file at %s' % os.environ['ARNIEFILE'])
        with open("%s" % os.environ["ARNIEFILE"], 'r') as f:
            for line in f.readlines():
                if line.strip():
                    if not line.startswith('#'):
                        key, string = line.split(':')
                        string = string.strip()
                        if key not in return_dct:
                            return_dct[key] = string

    # If the dict is empty, arnie couldn't find packages in the environment variables or arnie file
    if not return_dct:
        raise EnvironmentError("No prediction packages found in your environment. Check your environment variables or your ARNIEFILE.")

    # Some of the functions currently assume a TMP directory to save temporary files generated by various prediction packages.
    # If a TMP path is not provided, default to a folder in the current directory.
    if 'TMP' not in return_dct:
        tempPath = './tmp'
        if not os.path.exists(tempPath):
            os.mkdir(tempPath)
            
        return_dct['TMP'] = tempPath

    return return_dct


###############################################################################
# Structure representation conversion
###############################################################################

def convert_dotbracket_to_bp_list(s, allow_pseudoknots=False):
    if allow_pseudoknots:
        lefts = ["(", "[", "{", "<"]
        rights = [")", "]", "}", ">"]
        lower_alphabet = [chr(lower) for lower in range(97, 123)]
        upper_alphabet = [chr(upper) for upper in range(65, 91)]
        lefts.extend(lower_alphabet)
        rights.extend(upper_alphabet)
    else:
        lefts = ["("]
        rights = [")"]

    char_ignored = [char for char in s if char not in rights + lefts + ["."]]
    char_ignored = list(set(char_ignored))
    if char_ignored != []:
        print("WARNING: characters in structure,", char_ignored, "ignored!")

    l = []
    for left, right in zip(lefts, rights):
        bp1 = []
        bp2 = []
        for i, char in enumerate(s):
            if char == left:
                bp1.append(i)
            elif char == right:
                bp2.append(i)

        for i in list(reversed(bp1)):
            for j in bp2:
                if j > i:
                    l.append([i, j])

                    bp2.remove(j)
                    break
    l = sorted(l, key=lambda x: x[0])
    return l


def convert_dotbracket_to_bp_dict(s, allow_pseudoknots=False):

    if allow_pseudoknots:
        lefts = ["(", "[", "{", "<"]
        rights = [")", "]", "}", ">"]
        lower_alphabet = [chr(lower) for lower in range(97, 123)]
        upper_alphabet = [chr(upper) for upper in range(65, 91)]
        lefts.extend(lower_alphabet)
        rights.extend(upper_alphabet)
    else:
        lefts = ["("]
        rights = [")"]

    char_ignored = [char for char in s if char not in rights + lefts + ["."]]
    char_ignored = list(set(char_ignored))
    if char_ignored != []:
        print("WARNING: characters in structuture,", char_ignored, "ignored!")

    m = {}
    for left, right in zip(lefts, rights):
        bp1 = []
        bp2 = []
        for i, char in enumerate(s):
            if char == left:
                bp1.append(i)
            elif char == right:
                bp2.append(i)

        for i in list(reversed(bp1)):
            for j in bp2:
                if j > i:
                    m[i] = j
                    m[j] = i

                    bp2.remove(j)
                    break

    return m


def convert_dotbracket_to_matrix(s, allow_pseudoknots=False):
    matrix = np.zeros([len(s), len(s)])
    bp_list = convert_dotbracket_to_bp_dict(
        s, allow_pseudoknots=allow_pseudoknots)
    for k, v in bp_list.items():
        matrix[k, v] = 1
    return matrix


def convert_bp_list_to_dotbracket(bp_list, seq_len):
    db = "." * seq_len
    # group into bps that are not intertwined and can use same brackets!
    groups = _group_into_non_conflicting_bp(bp_list)

    # all bp that are not intertwined get (), but all others are
    # groups to be nonconflicting and then asigned (), [], {}, <> by group
    chars_set = [("(", ")"), ("(", ")"), ("[", "]"), ("{", "}"), ("<", ">")]
    alphabet = [(chr(lower), chr(upper))
                for upper, lower in zip(list(range(65, 91)), list(range(97, 123)))]
    chars_set.extend(alphabet)

    if len(groups) > len(chars_set):
        print("WARNING: PK too complex, not enough brackets to represent it.")

    for group, chars in zip(groups, chars_set):
        for bp in group:
            db = db[:bp[0]] + chars[0] + \
                db[bp[0] + 1:bp[1]] + chars[1] + db[bp[1] + 1:]
    return db


def get_bpp_from_dbn(dbn_struct):
    """ Gets a base-pairing matrix from a dot-bracket notation structure
    The matrix has 1's at position i,j if sequence positions i,j are base-paired
    and 0 otherwise.

    Args: 
      dbn_struct - structure in MFE format
    Return:
      2D numpy array
    """

    bpp_matrix = np.zeros((len(dbn_struct), len(dbn_struct)))

    bkt_open = "(<{"
    bkt_close = ")>}"

    for ii in range(len(bkt_open)):
        open_char = bkt_open[ii]
        close_char = bkt_close[ii]

        open_pos = []

        for jj in range(len(dbn_struct)):
            if dbn_struct[jj] == open_char:
                open_pos += [jj]
            if dbn_struct[jj] == close_char:
                pair = open_pos[-1]
                del open_pos[-1]
                bpp_matrix[jj, pair] = 1
                bpp_matrix[pair, jj] = 1

    return bpp_matrix

def get_helices(s, allowed_buldge_len=0):
    bp_list = convert_dotbracket_to_bp_list(s, allow_pseudoknots=True)
    helices = []
    current_helix = []
    while bp_list != []:
        current_bp = bp_list.pop(0)
        if current_helix == []:
            current_helix.append(current_bp)
        else:
            in_helix_left = list(range(current_helix[-1][0] + 1, current_helix[-1][0] + allowed_buldge_len + 2))
            in_helix_right = list(range(current_helix[-1][1] - allowed_buldge_len - 1, current_helix[-1][1]))
            if current_bp[0] in in_helix_left and current_bp[1] in in_helix_right:
                current_helix.append(current_bp)
            else:
                helices.append(current_helix)
                current_helix = [current_bp]
    helices.append(current_helix)
    return helices

def post_process_struct(s, allowed_buldge_len=0, min_len_helix=1):
    ''' given a structure, remove any helices that are too short

    allowed_buldge_len, if 0 does not allow any buldges when calculating the length
    of the helices, if 1 allows 0-1 and 1-1 buldges, if 2 additionally allows 2-0,2-1,2-2 bulges etc
    min_len_helix, any helices less than this value will be removed.
    '''
    helices = get_helices(s, allowed_buldge_len)
    bp_list_out = []
    for helix in helices:
        if len(helix) >= min_len_helix:
            bp_list_out.extend(helix)
    s_out = convert_bp_list_to_dotbracket(bp_list_out, len(s))
    return s_out


###############################################################################
# Evaluating a structure
###############################################################################

def get_expected_accuracy(dbn_string, bp_matrix, mode='mcc'):
    '''given a secondary structure as dbn string and base pair matrix, 
    assess expected accuracy for the structure.

    Inputs:
    dbn_string (str): Secondary structure string in dot-parens notation.
    bp_matrix (NxN array):  symmetric matrix of base pairing probabilities.
    mode: ['mcc','fscore','sen','ppv']: accuracy metric for which to compute expected value.

    Returns: expected accuracy value.
    '''

    assert bp_matrix.shape[0] == bp_matrix.shape[1]
    assert bp_matrix.shape[0] == len(dbn_string)

    struct_matrix = convert_dotbracket_to_matrix(dbn_string)
    N = len(dbn_string)

    pred_m = struct_matrix[np.triu_indices(N)]
    probs = bp_matrix[np.triu_indices(N)]

    TP = np.sum(np.multiply(pred_m, probs)) + 1e-6
    TN = 0.5*N*N-1 - np.sum(pred_m) - np.sum(probs) + TP + 1e-6
    FP = np.sum(np.multiply(pred_m, 1-probs)) + 1e-6
    FN = np.sum(np.multiply(1-pred_m, probs)) + 1e-6

    a, b = np.triu_indices(N)
    cFP = 1e-6  # compatible false positives
    # for i in range(len(pred_m)):
    #     if np.sum(struct_matrix,axis=0)[a[i]] + np.sum(struct_matrix,axis=0)[b[i]]==0:
    #        cFP += np.multiply(pred_m[i], 1-probs[i])

    if mode == 'sen':
        return TP/(TP + FN)
    elif mode == 'ppv':
        return TP/(TP + FP - cFP)
    elif mode == 'mcc':
        return (TP*TN - (FP - cFP)*FN)/np.sqrt((TP + FP - cFP)*(TP + FN)*(TN + FP - cFP)*(TN + FN))
    elif mode == 'fscore':
        return 2*TP/(2*TP + FP - cFP + FN)
    else:
        print('Error: mode not understood.')


def get_mean_base_pair_propensity(dbn_string):
    '''Measure of base pair locality.'''
    mat = convert_dotbracket_to_matrix(dbn_string)
    i, j = np.where(mat == 1)
    # divide by 2 because symmetric matrix
    mean_bp_dist = 0.5*np.mean([np.abs(x-y) for (x, y) in list(zip(i, j))])
    return mean_bp_dist

def is_PK(s):
    '''return if dotbracket structure represents a PK'''
    return ("[" in s) or ("{" in s) or ("<" in s)


def compare_structure_to_native(s, native, metric="all", PK_involved=None):
    if metric not in ["PPV", "sensitivity", "F1_score", "all"]:
        raise ValueError(
            'Only PPV, sensitivity, F1_score and all comparison currently implemented.')

    if PK_involved is None:
        s_list = convert_dotbracket_to_bp_list(s, allow_pseudoknots=True)
        native_list = convert_dotbracket_to_bp_list(
            native, allow_pseudoknots=True)
    elif PK_involved:
        s_list = _seperate_structure_into_PK_involved_or_not(s)["pk_bps"]
        native_list = _seperate_structure_into_PK_involved_or_not(native)[
            "pk_bps"]
    else:
        s_list = _seperate_structure_into_PK_involved_or_not(s)["no_pk_bps"]
        native_list = _seperate_structure_into_PK_involved_or_not(native)[
            "no_pk_bps"]

    PP = len(s_list)
    P = len(native_list)

    TP = len([x for x in s_list if x in native_list])  # true positives
    FP = len([x for x in s_list if x not in native_list])  # false positives
    FN = len([x for x in native_list if x not in s_list])  # false negative

    PPV = TP / PP if PP != 0 else 0
    sen = TP / P if P != 0 else 0
    F1 = (2 * PPV * sen) / (sen + PPV) if sen + PPV != 0 else 0

    if metric == "PPV":
        return PPV
    elif metric == "sensitivity":
        return sen
    elif metric == "F1_score":
        return F1
    elif metric == "all":
        return {"PPV": PPV, "sensitivity": sen, "F1_score": F1}


def compare_structures_to_natives(structs, natives, comparison="basepairs", metric="all"):
    # TODO other metrics? https://en.wikipedia.org/wiki/Positive_and_negative_predictive_values#Negative_predictive_value_(NPV)
    if comparison not in ["basepairs", "is_PK", "non_PK_basepairs", "PK_basepairs"]:
        raise ValueError(
            'Only basepairs and is_PK comparison currently implemented.')
    if metric not in ["PPV", "sensitivity", "F1_score", "all"]:
        raise ValueError(
            'Only PPV, sensitivity, F1_score and all comparison currently implemented.')

    if comparison == "basepairs":
        structs_list = [convert_dotbracket_to_bp_list(
            s, allow_pseudoknots=True) for s in structs]
        natives_list = [convert_dotbracket_to_bp_list(
            native, allow_pseudoknots=True) for native in natives]

        PP = sum(len(x) for x in structs_list)
        P = sum(len(x) for x in natives_list)

        TP = sum(len([z for z in x if z in y])
                 for x, y in zip(structs_list, natives_list))  # true positives

    elif comparison == "non_PK_basepairs":
        structs_list = [_seperate_structure_into_PK_involved_or_not(
            s)["no_pk_bps"] for s in structs]
        natives_list = [_seperate_structure_into_PK_involved_or_not(
            native)["no_pk_bps"] for native in natives]
        PP = sum(len(x) for x in structs_list)
        P = sum(len(x) for x in natives_list)

        TP = sum(len([z for z in x if z in y])
                 for x, y in zip(structs_list, natives_list))  # true positives

    elif comparison == "PK_basepairs":
        structs_list = [_seperate_structure_into_PK_involved_or_not(s)[
            "pk_bps"] for s in structs]
        natives_list = [_seperate_structure_into_PK_involved_or_not(
            native)["pk_bps"] for native in natives]
        PP = sum(len(x) for x in structs_list)
        P = sum(len(x) for x in natives_list)

        TP = sum(len([z for z in x if z in y])
                 for x, y in zip(structs_list, natives_list))  # true positives

    elif comparison == "is_PK":
        s_is_PK = [is_PK(s) for s in structs]
        native_is_PK = [is_PK(native) for native in natives]

        PP = sum(s_is_PK)
        P = sum(native_is_PK)

        # true positives
        TP = sum([x and y for x, y in zip(s_is_PK, native_is_PK)])
        # false positives
        FP = sum([x and not y for x, y in zip(s_is_PK, native_is_PK)])
        # false negatives
        FN = sum([not x and y for x, y in zip(s_is_PK, native_is_PK)])
        TN = sum([not x and not y for x, y in zip(
            s_is_PK, native_is_PK)])  # true negatives

    PPV = TP / PP if PP != 0 else 0
    sen = TP / P if P != 0 else 0
    F1 = (2 * PPV * sen) / (sen + PPV) if sen + PPV != 0 else 0

    if metric == "PPV":
        return PPV
    elif metric == "sensitivity":
        return sen
    elif metric == "F1_score":
        return F1
    elif metric == "all":
        return {"PPV": PPV, "sensitivity": sen, "F1_score": F1}


###############################################################################
# Structure helpers
###############################################################################


def _seperate_structure_into_PK_involved_or_not(s):
    bp_list = convert_dotbracket_to_bp_list(s, allow_pseudoknots=True)
    groups = _group_into_non_conflicting_bp(bp_list)
    bp_list_no_pk = groups[0]
    bp_list_pk = [bp for group in groups[1:] for bp in group]
    return {"no_pk_bps": bp_list_no_pk, "pk_bps": bp_list_pk}


def _get_non_redudant_bp_list(conflict_list):
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


def _get_list_bp_conflicts(bp_list):
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
                conflicts.append([current_bp, bp])
        return conflicts + _get_list_bp_conflicts(bp_list[1:])


def _group_into_non_conflicting_bp(bp_list):
    ''' given a bp_list, group basepairs into groups that do not conflict

    Args
            bp_list: list of base_pairs

    Returns:
            groups of baspairs that are not intertwined
    '''
    conflict_list = _get_list_bp_conflicts(bp_list)

    non_redudant_bp_list = _get_non_redudant_bp_list(conflict_list)
    bp_with_no_conflict = [
        bp for bp in bp_list if bp not in non_redudant_bp_list]
    groups = [bp_with_no_conflict]
    while non_redudant_bp_list != []:
        current_bp = non_redudant_bp_list[0]
        current_bp_conflicts = []
        for conflict in conflict_list:
            if current_bp == conflict[0]:
                current_bp_conflicts.append(conflict[1])
            elif current_bp == conflict[1]:
                current_bp_conflicts.append(conflict[0])
        max_group = [
            bp for bp in non_redudant_bp_list if bp not in current_bp_conflicts]
        to_remove = []
        for i, bpA in enumerate(max_group):
            for bpB in max_group[i:]:
                if bpA not in to_remove and bpB not in to_remove:
                    if [bpA, bpB] in conflict_list or [bpB, bpA] in conflict_list:
                        to_remove.append(bpB)
        group = [bp for bp in max_group if bp not in to_remove]
        groups.append(group)
        non_redudant_bp_list = current_bp_conflicts
        conflict_list = [conflict for conflict in conflict_list if conflict[0]
                         not in group and conflict[1] not in group]
    return groups


###############################################################################
# ORPHANED unused anywhere in package
# figure out utility, take outside of utils and document or depricate
###############################################################################


def convert_multiple_dbns_to_eternafold_input(seq, list_of_constraint_strings, filename):
    '''hard-coded to have 3 constraints right now for use in eternafold training with kd-ligand data.'''
    constraint_list = []
    for constraint_string in list_of_constraint_strings:
        constraint_list.append(write_constraint_string(seq, constraint_string))

    with open('%s' % filename, 'w') as out:
        for i in range(len(seq)):
            out.write('%d\t%s\t%d\t%d\t%d\n' % (
                i+1, seq[i], constraint_list[0][i], constraint_list[1][i], constraint_list[2][i]))


def write_reactivity_file_vienna(reactivities, sequence, fname=None):
    '''write reactivity (either SHAPE or DMS) to file format used by ViennaRNA.

      ex:
      1 U 0.120768
      2 G 0.190510
      3 U 0.155776

    Args:
      reactivities (list): a list of normalized reactivity float data. 
      sequence: RNA sequence
      Negative numbers can be used to indicate no signal.
    '''

    assert len(reactivities) == len(sequence)

    if fname is None:
        fname = '%s.SHAPE' % filename()

    with open(fname, 'w') as f:
        i = 1
        for char, reactivity in list(zip(sequence, reactivities)):
            if reactivity >= 0:
                f.write('%d %s %f\n' % (i, char, reactivity))

            i += 1
    return fname


def get_missing_motif_bases(seq):
    FMN_apt1 = 'AGGAUAU'
    FMN_apt2 = 'AGAAGG'
    a = seq.find(FMN_apt1) - 1
    # print(a)
    b = seq.find(FMN_apt2) + len(FMN_apt2)
    #print(seq.find(FMN_apt2), b)

    return seq[a], seq[b]


def combo_list_to_dbn_list(seq, final_combo_list, apt_idx_list, apt_ss_list):
    """ Helper function for write_constraints 
    Converts the combination of aptamer fragments to dbn strings
    """
    dbn_string_list = []
    for combo in final_combo_list:
        dbn_string = ['.']*len(seq)
        temp = apt_idx_list[combo, :][:, 1:3]

        # Fill in the dbn string with the aptamer ss
        for (start, finish), apt_ss in zip(temp, apt_ss_list):
            dbn_string[start:finish] = list(apt_ss)
        dbn_string = ''.join(dbn_string)

        # Check if aptamer is flipped
        if temp[0, 0] > temp[-1, 0]:
            dbn_string = flip_ss(dbn_string)
        dbn_string_list.append(dbn_string)
    return dbn_string_list


def write_combo_constraints(seq, raw_apt_seq, raw_apt_ss, verbose=False):
    """ Given a sequence, get all possible secondary constraints of the aptamer 

    Args:
      seq: RNA sequence
      raw_apt_seq: aptamer sequence
        e.g. CAAAG+CAAAG+GGCCUUUUGGCC
        + denotes splitable aptamer
      raw_apt_ss: aptamer secondary structure
        e.g. (xxx(+)xxx)+((((xxxx))))
        + denotes splitable aptamer
        x denotes unpaired base
        . denotes wildcard (can be anything)
      verbose: to be verbose
    Returns
      list of all possible dbn_string for the given aptamer
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
        raise ValueError(
            'Missing + in aptamer sequence or secondary structure')

    # Iterate through each aptamer fragment and save its idx,
    apt_idx_list = []
    for idx, (apt_seq, apt_ss) in enumerate(zip(apt_seq_list, apt_ss_list)):
        if seq.find(apt_seq) == -1:
            raise ValueError("Aptamer segment {} not found".format(idx+1))
        if len(apt_seq) != len(apt_ss):
            raise ValueError(
                "Mismatch between aptamer sequence and aptamer secondary structure")

        # Save all locations of each fragment
        for m in re.finditer(apt_seq, seq):
            # Note: cannot get overlapping segments
            span = m.span()
            start = span[0]
            finish = span[1]
            ss_length = len(apt_seq)
            apt_idx_list.append([idx, start, finish, ss_length])
    apt_idx_list = np.array(apt_idx_list)

    # Now combine aptamer fragments into full secondary structure constraints
    N_frag = len(apt_ss_list)  # Number of fragments to stitch together

    # Get a list of all possible combination (each combination is a list as well)
    temp = np.array([np.where(apt_idx_list[:, 0] == idx)[0]
                     for idx in range(N_frag)])
    if len(temp) > 1:
        combo = np.meshgrid(*temp)
        # reformat the combinations
        combo_list = np.array(combo).T.reshape(-1, N_frag)
    else:
        combo_list = np.array([[x] for x in temp[0]])

    # Check each combination to make sure it's feasible, if not remove
    final_combo_list = prune_combo_list(combo_list, apt_idx_list, N_frag)

    # Convert each combination into a dbn_string
    dbn_list = combo_list_to_dbn_list(
        seq, final_combo_list, apt_idx_list, apt_ss_list)

    # List will be empty if nothing can be in order (only for >= 3 fragments)
    return dbn_list


def write_constraint_string(seq, constraint_dbn):
    '''write set of integers to represent constraints, i.e. for use in bpseq format.'''

    assert(len(seq) == len(constraint_dbn))

    bp_list = convert_dotbracket_to_bp_dict(constraint_dbn)

    constraint_list = []

    for i, c in enumerate(constraint_dbn):
        if c == 'x':
            constraint = 0
        elif c == '.':
            constraint = -1  # or -1 if undefined
        elif c in ['(', ')']:
            constraint = bp_list[i]+1
        else:
            print('Error reading constraint string', c)
        constraint_list.append(constraint)
    return constraint_list


def prune_combo_list(combo_list, apt_idx_list, N_frag):
    """ Helper function for write_constraints 
    Prunes all possible combinations of the aptamer fragments
    """
    final_idx_list = []
    for idx, combo in enumerate(combo_list):
        temp = apt_idx_list[combo, :]
        start_idx_list = temp[:, 1]
        ss_len_list = temp[:, 3]
        if N_frag == 1:  # If <= 2 fragments than always in order
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
                if dx[0] > 0:  # fragments in increasing order
                    ss_len_list = ss_len_list[:-1]
                else:  # fragments in decreasing order
                    ss_len_list = ss_len_list[1:]

                # Check if there is enough bases to fit the aptamer fragments
                if np.all(np.abs(dx) >= ss_len_list):
                    final_idx_list.append(idx)
    final_combo_list = combo_list[final_idx_list]
    return np.array(final_combo_list)


def flip_ss(ss):
    """ Flips a secondary structure 
    Only flips the unpaired bases
    e.g. (.((....))....( ->
         ).((....))....)
    """
    bp_list = []
    unmatch_list = []
    ss = list(ss)
    for idx, c in enumerate(ss):
        if c == '(':
            bp_list.append(idx)
        elif c == ')':
            if len(bp_list):
                bp_list.pop()
            else:
                unmatch_list.append(idx)
    unmatch_list += bp_list
    for idx in unmatch_list:
        if ss[idx] == '(':
            ss[idx] = ')'
        else:
            ss[idx] = '('
    ss = ''.join(ss)
    return ss

def write_constraints(seq, motif=False, MS2=False, LIG=False, lig1=('nAGGAUAU', '(xxxxxx('), lig2=('AGAAGGn', ')xxxxx)')):
    '''Inputs:
    seq: RNA sequence
    motif: tuple (seq, struct) of motif. For example: PUM would be motif=('UGUAUAUA','xxxxxxxx').
    MS2: bool, whether to include MS2 constraint or not
    lig1: tuple (seq, struct) for 5' portion of ligand aptamer. Default is FMN.
    lig2: tuple (seq, struct) for 3' portion of ligand aptamer

    Outputs:
    dbn string, () for paired, x for unpaired, . for unspecified
    '''

    # when FMN aptamer and MS2 aptamer overlap, MS2 loses out on bp
    MS2_apt = 'ACAUGAGGAUCACCCAUGU'
    LIG_apt1 = lig1[0].replace('n', '')
    LIG_apt2 = lig2[0].replace('n', '')

    unpaired_list = []
    bp_list = {}

    dbn_string = ['.']*len(seq)

    if motif:
        if LIG:
            raise ValueError(
                'Sorry, due to some hacky hard-coding, cannot handle both motif and LIG inputs at this time.')
        else:
            return write_constraints(seq, LIG=True, lig1=motif, lig2=('', ''))

    if LIG:
        if seq.find(LIG_apt1) == -1:
            raise RuntimeError("ligand 5' aptamer domain not found, %s" % seq)

        else:
            if lig1[0].startswith('n'):
                start1 = seq.find(LIG_apt1) + len(LIG_apt1) - len(lig1[0])
                if start1 < 0:  # hws: changed from 1, maybe wrong?
                    start1 = seq.find(
                        LIG_apt1, start1+len(lig1[0])+1) + len(LIG_apt1) - len(lig1[0])
            else:
                start1 = seq.find(LIG_apt1)
                if start1 < 0:  # hws: changed from 1, maybe wrong?
                    start1 = seq.find(LIG_apt1, start1+len(lig1[0])+1)

            finish1 = start1 + len(lig1[0])

            if lig2[0].startswith('n'):
                start2 = seq.find(LIG_apt2, finish1+1) + \
                    len(LIG_apt2) - len(lig2[0])
            else:
                start2 = seq.find(LIG_apt2, finish1+1)
            finish2 = start2 + len(lig2[0])
            #print('start1, finish1, start2, finish2 FMN', start1, finish1, start2, finish2)
            dbn_string[start1:finish1] = list(lig1[1])
            dbn_string[start2:finish2] = list(lig2[1])

    if MS2:
        if seq.find(MS2_apt) == -1:
            raise RuntimeError("MS2 aptamer domain not found: %s" % seq)
        else:
            start = seq.find(MS2_apt)
            finish = start+len(MS2_apt)
            #print('start, finish MS2', start, finish)

            if dbn_string[start] != ".":
                #print('warning, aptamer overlap')
                dbn_string[start+1:finish-1] = list('((((x((xxxx))))))')
            else:
                dbn_string[start:finish] = list('(((((x((xxxx)))))))')

    return ''.join(dbn_string)