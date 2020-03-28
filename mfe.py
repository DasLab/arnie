import os, re, sys, shutil
import subprocess as sp
import random, string
import numpy as np
from .utils import *

DEBUG=False

# load package locations from yaml file, watch! global dict
package_locs = load_package_locations()

def mfe(seq, package='vienna_2', T=37,
    constraint=None, motif=None,
    dangles=True, noncanonical=False,
    bpps=False, param_file=None, coaxial=True, reweight=None,viterbi = False):
    ''' Compute MFE structure (within package) for RNA sequence.
    Note: this is distinct from the arnie MEA codebase, which takes any base pair probability matrix and computes the maximum expected accuracy structure.
    That said, Contrafold's default structure prediction is an MEA structure, not MFE.  In this module, calling Contrafold returns the default MEA structure unless the 
    --viterbi flag is used, which will do the viterbi (MFE) algorithm in contrafold.

        Args:
        seq (str): nucleic acid sequence
        T (float): temperature (Celsius)
        constraint (str): structure constraints
        motif (str): argument to vienna motif 
        dangles (bool): dangles or not, specifiable for vienna, nupack
        coaxial (bool): coaxial stacking or not, specifiable for rnastructure, vfold
        noncanonical(bool): include noncanonical pairs or not (for contrafold, RNAstructure (Cyclefold))

        Possible packages: 
        'vienna_2', 'vienna_1','contrafold_1','contrafold_2'
        
    Returns
        string: MFE structure
    '''

    try:
        pkg, version = package.lower().split('_')
    except:
        pkg, version = package.lower(), None

    if not bpps: # if bpps, already printed these warnings
        if not dangles and pkg not in ['vienna', 'nupack']:
            print('Warning: %s does not support dangles options' % pkg)
        if not coaxial and pkg not in ['rnastructure', 'vfold']:
            print('Warning: %s does not support coaxial options' % pkg)

    if pkg=='vienna':
        struct = mfe_vienna_(seq, version=version, T=T, dangles=dangles, constraint=constraint, motif=motif, param_file=param_file,reweight=reweight)
 
    elif pkg=='contrafold':
        struct = mfe_contrafold_(seq, version=version, T=T, constraint=constraint, param_file=param_file,viterbi=viterbi)

    else:
        raise ValueError('package %s not understood.' % package)

    return struct

def mfe_vienna_(seq, T=37, version='2', constraint=None, motif=None, param_file=None, dangles=True, reweight=None):
    """get partition function structure representation and Z

    Args:
        seq (str): nucleic acid sequence
        T (float): temperature
        constraint (str): structure constraints
        motif (str): argument to vienna motif  
    Returns
        str, float: secondary structure representation and Z
    """

    if not version:
        version='2'

    if version.startswith('2'):
        LOC=package_locs['vienna_2']
    elif version.startswith('1'):
        LOC=package_locs['vienna_1']
    else:
        raise RuntimeError('Error, vienna version %s not present' % version)

    command = ['%s/RNAfold' % LOC, '-T', str(T), '-p0'] #p0 doesn't predict bpps, saves time
    if motif is not None:
        command.append('--motif="%s"' % motif)

    if constraint is not None:
        fname = write([seq, constraint])
        command.append('-C')
        command.append('--enforceConstraint')
    else:
        fname = write([seq])

    if not dangles:
        command.append('--dangles=0')
        
    if reweight is not None:
        command.append('--commands=%s' % reweight)

    if param_file:
        command.append('--paramFile=%s' % param_file)

    with open(fname) as f:
        if DEBUG: print(fname)
        if DEBUG: print(' '.join(command))
        p = sp.Popen(command, stdin=f, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()

    if DEBUG:
        print('stdout')
        print(stdout)
        print('stderr')
        print(stderr)

    if p.returncode:
        raise Exception('RNAfold failed: on %s\n%s' % (seq, stderr))
    os.remove(fname)
    os.remove('rna.ps')

    if 'omitting constraint' in stderr.decode('utf-8'):
        raise ValueError('Constraint caused impossible structure')
    else:
        return stdout.decode('utf-8').split('\n')[1].split(' ')[0]


def mfe_contrafold_(seq, T=37, version='2', constraint=None, param_file=None,viterbi=False):
    """get partition function structure representation and free energy

    Args:
        seq (str): nucleic acid sequence
        T (float): temperature
        constraint (str): structure constraints
        motif (str): argument to vienna motif  
    Returns
        float: partition function
        Note: If the constraint is impossible then Z wil be equal to the Z unconstrained
    """
    if not version: version='2'

    fname = '%s.in' % filename()

    if version.startswith('2'):
        LOC=package_locs['contrafold_2']
    elif version.startswith('1'):
        LOC=package_locs['contrafold_1']
    else:
        raise RuntimeError('Error, Contrafold version %s not present' % version)

    command = ['%s/contrafold' % LOC, 'predict', fname]

    if param_file is not None:
        command = command + ['--params', param_file]

    if viterbi:
        command.append('--viterbi')

    if constraint is not None:
        convert_dbn_to_contrafold_input(seq, constraint, fname)
        command.append('--constraints')
    else:
        convert_dbn_to_contrafold_input(seq, ''.join(['.' for x in range(len(seq))]), fname)

    if DEBUG: print(' '.join(command))

    p = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE)

    stdout, stderr = p.communicate()

    if DEBUG:
        print('stdout')
        print(stdout)
        print('stderr')
        print(stderr)
    if p.returncode:
        raise Exception('Contrafold failed: on %s\n%s' % (seq, stderr))

    os.remove(fname)
    
    return stdout.decode('utf-8').split('\n')[-2]
