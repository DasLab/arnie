import os, re, sys, shutil
import subprocess as sp
import random, string
import numpy as np
from .utils import *

DEBUG=False

# load package locations from yaml file, watch! global dict
package_locs = load_package_locations()

def sample_structures(seq, n_samples = 10, package='vienna_2', T=37, constraint=None, 
	dangles=True, reweight=None, nonredundant=False):
    ''' Draw stochastic sampled structures for RNA sequence. Possible packages: 'vienna_1', 'vienna_2'

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

    else:
        raise ValueError('package %s either not understood or not supported at this moment.' % package)

    return struct_list

def sample_vienna_(seq, n_samples=10, T=37, version='2', constraint=None, 
	dangles=True, reweight=None, nonredundant=False):
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

    command = ['%s/RNAsubopt' % LOC, '-T', str(T), '--stochBT_en=%d' % n_samples]#, '-N']

    if constraint is not None:
        fname = write([seq, constraint])
        command.append('-C')
        command.append('--enforceConstraint')
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
            struct_list.append(line.split(' ')[0].replace('.','x')) #explicitly converting .'s to x's here to maintain x=unpaired, .=unconstrained
            # prob_list.append(float(line.split(' ')[-2]))
            # energy_list.append(float(line.split(' ')[-1]))

    return struct_list
