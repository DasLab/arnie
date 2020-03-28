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
    out.write('SS:\n')
    for x in SS_list:
      out.write('%d\n' % x)
    out.write('-1\nPairs:\n')
    for x,y in pairs_list:
      out.write('%d %d\n' % (x,y))
    out.write('-1 -1')

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

def print_available_packages():
  package_dct = load_package_locations()
  for key,v in package_dct.items():
    if key != "TMP":
      print(key,v)

def package_list():
  pkg_list=[]
  package_dct = load_package_locations()
  for key,v in package_dct.items():
    if key != "TMP":
      if v != "None":
        pkg_list.append(key)
  return pkg_list

def load_package_locations():
    '''Read in user-supplied file to specify paths to RNA folding packages. Specify this in your ~/.bashrc as $ARNIEFILE'''
    return_dct={}
    package_path = os.path.dirname(arnie.__file__)
    with open("%s" % os.environ["ARNIEFILE"],'r') as f:
        for line in f.readlines():
            if line.strip():
              if not line.startswith('#'):
                key, string = line.split(':')
                string = string.strip()
                return_dct[key] = string
    return return_dct
