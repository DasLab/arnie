import os, re, sys, shutil
import subprocess as sp
import random, string
import numpy as np
from .utils import *

DEBUG=False

# load package locations from yaml file, watch! global dict
package_locs = load_package_locations()

def sample_structures(seq, n_samples = 10, package='vienna_2', T=37, constraint=None, param_file=None,
	dangles=True, reweight=None, nonredundant=False):
    ''' Draw stochastic sampled structures for RNA sequence. Possible packages: 'eternafold', 'vienna_2'

        Args:
        seq (str): nucleic acid sequence
        T (float): temperature (Celsius)
        constraint (str): structure constraints
        motif (str): argument to vienna motif 
        dangles (bool): dangles or not, specifiable for vienna, nupack
        noncanonical(bool): include noncanonical pairs or not (for contrafold, RNAstructure (Cyclefold))
        
    Returns
        list of structures
        list of energies
        list of probabilities 
    '''

    try:
        pkg, version = package.lower().split('_')
    except:
        pkg, version = package.lower(), None

    if not dangles and pkg not in ['vienna','nupack']:
        print('Warning: %s does not support dangles options' % pkg)

    if pkg=='vienna':
        struct_list = sample_vienna_(seq, n_samples=n_samples, version=version, T=T, 
        	dangles=dangles, constraint=constraint, reweight=reweight, nonredundant = nonredundant)

    elif pkg=='eternafold':
        struct_list = sample_eternafold_(seq, n_samples=n_samples, param_file=param_file, constraint=constraint, nonredundant = nonredundant)

    else:
        raise ValueError('package %s either not understood or not supported at this moment.' % package)

    return struct_list

def sample_vienna_(seq, n_samples=10, T=37, version='2', constraint=None, 
	dangles=True, reweight=None, nonredundant=False):
    """Stochastically sample structures from Vienna RNAsubopt.

    Inputs:
        seq (str): nucleic acid sequence
        n_samples (int): number of structures to sample.
        T (float): temperature
        constraint (str): structure constraints
        motif (str): argument to vienna motif  
    Outputs:
        struct_list (list): list of stochastically-sampled structures.
    """

    if not version:
        version='2'

    if version.startswith('2'):
        LOC=package_locs['vienna_2']
    elif version.startswith('1'):
        LOC=package_locs['vienna_1']
    else:
        raise RuntimeError('Error, vienna version %s not present' % version)

    command = ['%s/RNAsubopt' % LOC, '-T', str(T), '--stochBT_en=%d' % n_samples]#, '-N']

    if constraint is not None:
        fname = write([seq, constraint])
        command.append('-C')
        #command.append('--enforceConstraint')
    else:
        fname = write([seq])

    if not dangles:
        command.append('--dangles=0')

    if nonredundant:
    	command.append('-N')

    if reweight is not None:
        command.append('--commands=%s' % reweight)

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
        raise Exception('RNAsubopt failed: on %s\n%s' % (seq, stderr))
    os.remove(fname)

    if 'omitting constraint' in stderr.decode('utf-8'):
        raise RuntimeError("Constraint omitted, Impossible structure")

    else:
        struct_list, prob_list, energy_list = [],[],[]
        output_lines = stdout.decode('utf-8').split('\n')[1:-1] # first line is just repeating sequence, last is empty space
        for line in output_lines:
            struct_list.append(line.split(' ')[0])
            # prob_list.append(float(line.split(' ')[-2]))
            # energy_list.append(float(line.split(' ')[-1]))

    return struct_list

def sample_eternafold_(seq, n_samples=10, param_file=None, constraint=None, nonredundant=False):
    """Stochastically sample structures from EternaFold.

    Inputs:
        seq (str): nucleic acid sequence
        n_samples (int): number of structures to sample.
        T (float): temperature
        constraint (str): structure constraints
        motif (str): argument to vienna motif  
    Outputs:
        struct_list (list): list of stochastically-sampled structures.
    """

    fname = '%s.in' % filename()
    LOC=package_locs['eternafold']


    command = ['%s/contrafold' % LOC, 'sample', fname]

    if param_file is not None:
        command = command + ['--params', param_file]
    else:
        command = command + ['--params', package_locs['eternafoldparams']]

    if constraint is not None:
        convert_dbn_to_contrafold_input(seq, constraint, fname)
        command.append('--constraints')
    else:
        convert_dbn_to_contrafold_input(seq, ''.join(['.' for x in range(len(seq))]), fname)

    if DEBUG: print(' '.join(command))

    p = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE)

    stdout, stderr = p.communicate()

    struct_list = stdout.decode('utf-8').split('\n')[:-1]

    if DEBUG:
        print('stdout')
        print(stdout)
        print('stderr')
        print(stderr)
    if p.returncode:
        raise Exception('Eternafold sample failed: on %s\n%s' % (seq, stderr))

    os.remove(fname)
    return struct_list
