import os, re
import subprocess as sp
import random, string
import numpy as np
import arnie

def write_vector_to_file(vector, outfile):
  for x in vector:
    outfile.write('%.3f\n' % x)
  return

def write_matrix_to_file(matrix, outfile):
  for x in matrix:
    outfile.write('\t'.join(['%.3f' % y for y in x])+'\n')
  return

def complement_to_(string):
        base_pairing_dct = {'a':'u', 'u':'a', 'g':'c', 'c':'g','t':'a'}
        return ''.join(base_pairing_dct[x.lower()] for x in string[::-1])

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

  a,b = np.triu_indices(N)
  cFP = 1e-6 # compatible false positives
  # for i in range(len(pred_m)):
  #     if np.sum(struct_matrix,axis=0)[a[i]] + np.sum(struct_matrix,axis=0)[b[i]]==0:
  #        cFP += np.multiply(pred_m[i], 1-probs[i])

  if mode=='sen':
    return TP/(TP + FN)
  elif mode=='ppv':
    return TP/(TP + FP - cFP)
  elif mode=='mcc':
    return (TP*TN - (FP - cFP)*FN)/np.sqrt((TP + FP - cFP)*(TP + FN)*(TN + FP - cFP)*(TN + FN))
  elif mode=='fscore':
    return 2*TP/(2*TP + FP - cFP + FN)
  else:
    print('Error: mode not understood.')


def get_mean_base_pair_propensity(dbn_string):
    '''Measure of base pair locality.'''
    mat = convert_dotbracket_to_matrix(dbn_string)
    i, j = np.where(mat==1)
    #divide by 2 because symmetric matrix
    mean_bp_dist = 0.5*np.mean([np.abs(x-y) for (x,y) in list(zip(i,j))])
    return mean_bp_dist

def convert_dotbracket_to_bp_list(s, allow_pseudoknots=False):
    m = {}
    bp1=[]
    bp2=[]
    bp1_pk=[]
    bp2_pk=[]
    for i, char in enumerate(s):
        if char=='(':
            bp1.append(i)
        if char==')':
            bp2.append(i)
        if allow_pseudoknots:
            if char=='[':
              bp1_pk.append(i)
            if char==']':
              bp2_pk.append(i)

    for i in list(reversed(bp1)):
        for j in bp2:
            if j > i:
                m[i]=j
                m[j]=i

                bp2.remove(j)
                break
    for i in list(reversed(bp1_pk)):
        for j in bp2_pk:
            if j > i:
                m[i]=j
                m[j]=i

                bp2_pk.remove(j)
                break
    return m

def convert_dotbracket_to_matrix(s, allow_pseudoknots=False):
  matrix=np.zeros([len(s),len(s)])
  bp_list = convert_dotbracket_to_bp_list(s, allow_pseudoknots=allow_pseudoknots)
  for k,v in bp_list.items():
    matrix[k,v] = 1
  return matrix


def convert_dbn_to_RNAstructure_input(seq, constraints, filename):
  assert(len(seq) == len(constraints))

  bp_list = convert_dotbracket_to_bp_list(constraints)

  SS_list, pairs_list = [], []

  for i, (s, c) in enumerate(list(zip(seq, constraints))):
    if c=='x':
      SS_list.append(i+1)
    elif c=='.':
      pass
    elif c =='(':
      pairs_list.append([i+1, bp_list[i]+1])
    elif c ==')':
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
    for x,y in pairs_list:
      out.write('%d %d\n' % (x,y))
    out.write('-1 -1')
    out.write('FMN:\n-1\nForbids:\n-1\n')

def write_constraint_string(seq, constraint_dbn):
  '''write set of integers to represent constraints, i.e. for use in bpseq format.'''

  assert(len(seq) == len(constraint_dbn))

  bp_list = convert_dotbracket_to_bp_list(constraint_dbn)

  constraint_list = []

  for i, c in enumerate(constraint_dbn):
        if c=='x':
          constraint=0
        elif c=='.':
          constraint=-1 #or -1 if undefined
        elif c in ['(',')']:
          constraint=bp_list[i]+1
        else:
          print('Error reading constraint string', c)
        constraint_list.append(constraint)
  return constraint_list

def convert_dbn_to_contrafold_input(seq, constraints, filename):
  constraint_list = write_constraint_string(seq, constraints)
  with open('%s' % filename, 'w') as out:
    for i in range(len(seq)):
      out.write('%d\t%s\t%d\n'%(i+1, seq[i], constraint_list[i]))

def convert_multiple_dbns_to_eternafold_input(seq, list_of_constraint_strings, filename):
  '''hard-coded to have 3 constraints right now for use in eternafold training with kd-ligand data.'''
  constraint_list=[]
  for constraint_string in list_of_constraint_strings:
    constraint_list.append(write_constraint_string(seq, constraint_string))
    
  with open('%s' % filename, 'w') as out:
    for i in range(len(seq)):
      out.write('%d\t%s\t%d\t%d\t%d\n' % (i+1, seq[i], constraint_list[0][i], constraint_list[1][i], constraint_list[2][i]))

def write_constraints(seq, motif=False, MS2=False, LIG=False, lig1=('nAGGAUAU','(xxxxxx('), lig2=('AGAAGGn',')xxxxx)')):
  '''Inputs:
  seq: RNA sequence
  motif: tuple (seq, struct) of motif. For example: PUM would be motif=('UGUAUAUA','xxxxxxxx').
  MS2: bool, whether to include MS2 constraint or not
  lig1: tuple (seq, struct) for 5' portion of ligand aptamer. Default is FMN.
  lig2: tuple (seq, struct) for 3' portion of ligand aptamer

  Outputs:
  dbn string, () for paired, x for unpaired, . for unspecified
  '''

  #when FMN aptamer and MS2 aptamer overlap, MS2 loses out on bp
  MS2_apt='ACAUGAGGAUCACCCAUGU'
  LIG_apt1=lig1[0].replace('n','')
  LIG_apt2=lig2[0].replace('n','')

  unpaired_list=[]
  bp_list={}

  dbn_string=['.']*len(seq)

  if motif:
    if LIG:
      raise ValueError('Sorry, due to some hacky hard-coding, cannot handle both motif and LIG inputs at this time.')
    else:
      return write_constraints(seq,LIG=True, lig1=motif,lig2=('',''))

  if LIG:
      if seq.find(LIG_apt1) == -1:
        raise RuntimeError("ligand 5' aptamer domain not found, %s" % seq)

      else:
        if lig1[0].startswith('n'):
          start1 = seq.find(LIG_apt1) + len(LIG_apt1) - len(lig1[0])
          if start1 < 0: #hws: changed from 1, maybe wrong?
            start1 = seq.find(LIG_apt1,start1+len(lig1[0])+1) + len(LIG_apt1) - len(lig1[0])
        else:
          start1 = seq.find(LIG_apt1)
          if start1 < 0: #hws: changed from 1, maybe wrong?
            start1 = seq.find(LIG_apt1,start1+len(lig1[0])+1)

        finish1 = start1 + len(lig1[0])   

        if lig2[0].startswith('n'):       
          start2 = seq.find(LIG_apt2, finish1+1) + len(LIG_apt2) - len(lig2[0])
        else:
          start2 = seq.find(LIG_apt2, finish1+1)
        finish2 = start2 + len(lig2[0])         
        #print('start1, finish1, start2, finish2 FMN', start1, finish1, start2, finish2)
        dbn_string[start1:finish1] = list(lig1[1])
        dbn_string[start2:finish2] = list(lig2[1])

  if MS2:
      if seq.find(MS2_apt) == -1:
        raise RuntimeError("MS2 aptamer domain not found: %s" %seq)
      else:
        start=seq.find(MS2_apt)
        finish=start+len(MS2_apt)
        #print('start, finish MS2', start, finish)

        if dbn_string[start] != ".":
          #print('warning, aptamer overlap')
          dbn_string[start+1:finish-1]=list('((((x((xxxx))))))')
        else:
          dbn_string[start:finish]=list('(((((x((xxxx)))))))')


  return ''.join(dbn_string)

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
        bpp_matrix[jj,pair] = 1
        bpp_matrix[pair,jj] = 1

  return bpp_matrix

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
    
    # Check if aptamer is flipped
    if temp[0,0] > temp[-1,0]:
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
    raise ValueError('Missing + in aptamer sequence or secondary structure')

  # Iterate through each aptamer fragment and save its idx,
  apt_idx_list = []
  for idx, (apt_seq, apt_ss) in enumerate(zip(apt_seq_list, apt_ss_list)):
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

def get_missing_motif_bases(seq):
  FMN_apt1='AGGAUAU'
  FMN_apt2='AGAAGG'
  a = seq.find(FMN_apt1) - 1
  #print(a)
  b = seq.find(FMN_apt2) + len(FMN_apt2)
  #print(seq.find(FMN_apt2), b)

  return seq[a], seq[b]

def run_RNAPVmin(probing_signal, seq, LOC, DEBUG, tauSigmaRatio=1, shapeConversion='S'):
    reac_file = write_reactivity_file_vienna(probing_signal, seq)
    fname = write([seq])
    RNApvmin_command = ['%s/RNApvmin' % LOC, reac_file, '--shapeConversion=%s' % shapeConversion, '--tauSigmaRatio=%f' % tauSigmaRatio]

    with open(fname) as f:
        if DEBUG: print(fname)
        if DEBUG: print(' '.join(RNApvmin_command))
        p = sp.Popen(RNApvmin_command, stdin=f, stdout=sp.PIPE, stderr=sp.PIPE)
    rnapvmin_stdout, rnapvmin_stderr = p.communicate()

    shape_file = filename()

    with open(shape_file,'wb') as f:
        f.write(rnapvmin_stdout)

    if DEBUG:
        print('stdout')
        print(rnapvmin_stdout)
        print('stderr')
        print(rnapvmin_stderr)

    if p.returncode:
        raise Exception('RNApvmin failed: on %s\n%s' % (seq, stderr))

    return shape_file

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

  with open (fname, 'w') as f:
    i = 1
    for char, reactivity in list(zip(sequence, reactivities)):
      if reactivity >= 0:
        f.write('%d %s %f\n' % (i, char, reactivity))

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

  with open (fname, 'w') as f:
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

def print_path_files():
  package_dct = load_package_locations()
  for key,v in package_dct.items():
    print(key,v)

def package_list():
  pkg_list=[]
  package_dct = load_package_locations()
  for key,v in package_dct.items():
    if key != "TMP" and key.lower() != 'bprna':
      if not key.startswith('linear'):
        if key == 'eternafoldparams' and 'eternafold' not in pkg_list:
          pkg_list.append('eternafold')
        else:
          if v != "None":
            pkg_list.append(key)
  return pkg_list

def load_package_locations(DEBUG=False):
    '''Read in user-supplied file to specify paths to RNA folding packages. Specify this in your ~/.bashrc as $ARNIEFILE'''
    return_dct={}
    package_path = os.path.dirname(arnie.__file__)

    if DEBUG: print('Reading Arnie file at %s' % os.environ['ARNIEFILE'])
    
    with open("%s" % os.environ["ARNIEFILE"],'r') as f:
      for line in f.readlines():
          if line.strip():
            if not line.startswith('#'):
              key, string = line.split(':')
              string = string.strip()
              return_dct[key] = string

    # if 'eternafoldparams' not in return_dct.keys():
    #   if 'eternafold' in return_dct.keys():                                         
    #     return_dct['eternafoldparams'] = "%s/../parameters/EternaFoldParams.v1" % return_dct['eternafold']

    return return_dct
