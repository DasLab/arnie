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
    linear=False, return_dG_MFE = False,
    dangles=True, noncanonical=False, beam_size=100,
    bpps=False, param_file=None, coaxial=True, reweight=None,viterbi = False,
    probing_signal=None,probing_kws=None, pseudo=False,
    shape_signal=None, dms_signal=None, shape_file=None, dms_file=None, **kwargs):

    ''' Compute MFE structure (within package) for RNA sequence.
    Note: this is distinct from the arnie MEA codebase, which takes any base pair probability matrix and computes the maximum expected accuracy structure.
    That said, Contrafold's default structure prediction is an MEA structure, not MFE.  In this module, calling Contrafold returns the default MEA structure unless the 
    --viterbi flag is used, which will do the viterbi (MFE) algorithm in contrafold.

        Args:
        seq (str): nucleic acid sequence
        T (float): temperature (Celsius)
        constraint (str): structure constraints
        linear (bool): call LinearFold to estimate MFE in Vienna or Contrafold
        motif (str): argument to vienna motif 
        return_dG_MFE (bool): also return dG(MFE) (specific to linearfold)
        dangles (bool): dangles or not, specifiable for vienna, nupack
        coaxial (bool): coaxial stacking or not, specifiable for rnastructure, vfold
        noncanonical(bool): include noncanonical pairs or not (for contrafold, RNAstructure (Cyclefold))
        shape_signal(list): list of normalized SHAPE reactivities, with negative values indicating no signal
        dms_signal(list): list of normalized DMS reactivities, with negative values indicating no signal
        pseudo: if True, will predict pseudoknots

        Possible packages: 
        'vienna_2', 'vienna_1','contrafold_1','contrafold_2', 'rnastructure'
        
    Returns
        string: MFE structure
    '''

    # TODO: update to be just probing_signal
    # if shape_signal:
    #     print('Warning: shape_signal is deprecated, use probing_signal')
    #     probing_signal = copy(shape_signal)

    # if dms_signal:
    #     print('Warning: dms_signal is deprecated, use probing_signal')
    #     probing_signal = copy(dms_signal)

    try:
        pkg, version = package.lower().split('_')
    except:
        pkg, version = package.lower(), None

    if not bpps: # if bpps, already printed these warnings
        if not dangles and pkg not in ['vienna', 'nupack']:
            print('Warning: %s does not support dangles options' % pkg)
        if not coaxial and pkg not in ['rnastructure', 'vfold']:
            print('Warning: %s does not support coaxial options' % pkg)

    if linear and pkg not in ['vienna','contrafold','eternafold']:
        print('Warning: LinearFold only implemented for vienna, contrafold, eternafold.')

    if pseudo and not pkg in ['rnastructure', 'nupack']:
        print('Warning: %s and pseudoknots not supported in Arnie yet' % pkg)

    if pkg=='vienna':
        if linear:
            if return_dG_MFE:
                struct, dG_MFE = mfe_linearfold_(seq, package='vienna', return_dG_MFE=return_dG_MFE, beam_size=beam_size)
            else:
                struct = mfe_linearfold_(seq, package='vienna', return_dG_MFE=return_dG_MFE)
        else:
            struct = mfe_vienna_(seq, version=version, T=T, dangles=dangles, constraint=constraint, motif=motif, param_file=param_file,
                reweight=reweight, probing_signal=probing_signal, **kwargs)
 
    elif pkg=='contrafold':
        if linear:
            if return_dG_MFE:
                struct, dG_MFE = mfe_linearfold_(seq, package='contrafold', return_dG_MFE=return_dG_MFE, beam_size=beam_size)
            else:
                struct = mfe_linearfold_(seq, package='contrafold', return_dG_MFE=return_dG_MFE)
        else:
            struct = mfe_contrafold_(seq, version=version, T=T, constraint=constraint, param_file=param_file,viterbi=viterbi)

    elif pkg=='eternafold':
        if linear:
            if return_dG_MFE:
                struct, dG_MFE = mfe_linearfold_(seq, package='eternafold', return_dG_MFE=return_dG_MFE, beam_size=beam_size)
            else:
                struct = mfe_linearfold_(seq, package='eternafold', return_dG_MFE=return_dG_MFE)

        else:

            if 'eternafoldparams' in package_locs.keys() and 'eternafold' not in package_locs.keys():
                struct = mfe_contrafold_(seq, version=version, T=T, constraint=constraint, param_file=package_locs['eternafoldparams'],viterbi=viterbi)

            elif 'eternafold' in package_locs.keys():
                
                # Using eternafold code and params in eternafold codebase
                efold_param_file = os.environ['ETERNAFOLD_PARAMETERS'] if os.environ.get('ETERNAFOLD_PARAMETERS') else package_locs['eternafold']+'/../parameters/EternaFoldParams.v1'
                if not os.path.exists(efold_param_file):
                    raise RuntimeError('Error: Parameters not found at %s' % efold_param_file)
                else:
                    struct = mfe_contrafold_(seq, version=version, T=T, constraint=constraint, DIRLOC=package_locs['eternafold'],
                        param_file=efold_param_file,viterbi=viterbi, probing_signal=probing_signal, probing_kws=probing_kws)


    elif pkg=='rnastructure':
        if linear:
            raise ValueError('package %s is not supported with linearfold.' % package)
        else:
            struct = mfe_rnastructure_(seq, version=version, T=T, constraint=constraint, 
                probing_signal=probing_signal,
                param_file=param_file, shape_signal=shape_signal, dms_signal=dms_signal, 
                shape_file=shape_file, dms_file=dms_file, pseudo = pseudo)
    else:
        raise ValueError('package %s not understood.' % package)

    if return_dG_MFE:
        return struct, dG_MFE
    else:
        return struct

def mfe_vienna_(seq, T=37, version='2', constraint=None, motif=None, param_file=None, dangles=True, reweight=None,
    probing_signal=None, shapeMethod='W', probing_kws=None, **kwargs):
    """get minimum free energy structure with Vienna

    Args:
        seq (str): nucleic acid sequence
        T (float): temperature
        constraint (str): structure constraints
        motif (str): argument to vienna motif  
    Returns
        str: secondary structure representation for MFE
    """

    if not version:
        version='2'

    if version.startswith('2'):
        LOC=package_locs['vienna_2']
    elif version.startswith('1'):

        LOC=package_locs['vienna_1']

    else:
        raise RuntimeError('Error, vienna version %s not present' % version)

    if constraint is not None:
        fname = write([seq, constraint])
        command.append('-C')
        command.append('--enforceConstraint')
    else:
        fname = write([seq])

    command = ['%s/RNAfold' % LOC, '-T', str(T), '-p0'] #p0 doesn't predict bpps, saves time
    if motif is not None:
        command.append('--motif=%s' % motif)

    if probing_signal is not None:
        if probing_kws is None:
            probing_kws={}

        if shapeMethod=='W':
            probing_file = run_RNAPVmin(probing_signal, seq, LOC, DEBUG, **probing_kws)

        elif shapeMethod=='D' or shapeMethod=='Z':
            probing_file = write_reactivity_file_RNAstructure(probing_signal)
            command.append('--shapeConversion=O')

        command.append('--shape=%s' % probing_file)
        command.append('--shapeMethod=%s' % shapeMethod)

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
    if os.path.exists('rna.ps'):
        os.remove('rna.ps')

    if 'omitting constraint' in stderr.decode('utf-8'):
        raise ValueError('Constraint caused impossible structure')
    else:
        return stdout.decode('utf-8').split('\n')[1].split(' ')[0]

def mfe_rnastructure_(seq, T=24, version=None, constraint=None, param_file=None, probing_signal=None,probing_kws=None,
    shape_signal=None, dms_signal=None, shape_file=None, dms_file=None, pseudo=False):
    """get minimum free energy structure
        with SHAPE or DMS data, uses the default slope and intercept in RNAStructure

    Args:
        seq (str): nucleic acid sequence
        T (float): temperature
    Returns
        float: MFE structure
    """

    if probing_signal is not None:
        shape_signal = probing_signal

    if param_file is not None:
        raise ValueError('Cannot run RNAstructure with non-default RNA parameters as specified in: %s' % param_file)
    if version is not None:
        raise ValueError('Cannot run RNAstructure with non-default version: %s' % version)
    if (shape_signal is not None) and (shape_file is not None):
        raise ValueError('Please specify SHAPE reactivities either as a list or in a SHAPE reactivity file')
    if (dms_signal is not None) and (dms_file is not None):
        raise ValueError('Please specify DMS reactivities either as a list or in a DMS reactivity file')

    LOC=package_locs['rnastructure']

    seq_file = write(['>sequence', seq])
    ct_fname = '%s.ct' % filename()

    command = []
    if not pseudo:
        command = command + ['%s/Fold' % LOC, seq_file, ct_fname, '-T', str(T + 273.15)]
    else:
        command = command + ['%s/ShapeKnots' % LOC, seq_file, ct_fname]
        # if dms_signal is not None:
        #     raise ValueError('Cannot run RNAstructure with DMS signal and pseudoknots.')
        if constraint is not None:
            raise ValueError('Cannot run RNAstructure with constraints and pseudoknots.')
    
    con_fname = None
    dms_fname = None
    shape_fname = None

    if constraint is not None:
        con_fname = '%s.CON' % filename()
        convert_dbn_to_RNAstructure_input(seq, constraint, con_fname)
        command.extend(['--constraint', con_fname])

    if dms_signal is not None:
        if len(dms_signal) != len(seq):
            raise RuntimeError('DMS signal used with RNAstructure must have same length as the sequence.')
        dms_fname = write_reactivity_file_RNAstructure(dms_signal)
        command.extend(['--DMS', dms_fname])

    if dms_file is not None:
        command.extend(['--DMS', dms_file])

    if shape_signal is not None:
        if len(shape_signal) != len(seq):
            raise RuntimeError('SHAPE signal used with RNAstructure must have same length as the sequence.')
        shape_fname = write_reactivity_file_RNAstructure(shape_signal)
        command.extend(['--SHAPE', shape_fname])

    if shape_file is not None:
        command.extend(['--SHAPE', shape_file])

    if DEBUG: print(' '.join(command))

    p = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE)

    stdout, stderr = p.communicate()

    if DEBUG:
        print('stdout')
        print(stdout)
        print('stderr')
        print(stderr)
    if p.returncode:
        raise Exception('RNAstructure failed: on %s\n%s' % (seq, stderr))

    if con_fname is not None:
        os.remove(con_fname)
    if dms_fname is not None:
        os.remove(dms_fname)
    if shape_fname is not None:
        os.remove(shape_fname)
    if seq_file is not None:
        os.remove(seq_file)

    dot_fname = '%s.dbn' % filename()
    command = ['%s/ct2dot' % LOC, ct_fname, "1", dot_fname]

    if DEBUG: print(' '.join(command))

    p = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE)

    stdout, stderr = p.communicate()

    if DEBUG:
        print('stdout')
        print(stdout)
        print('stderr')
        print(stderr)
    if p.returncode:
        raise Exception('RNAstructure ct2dot failed: on %s\n%s' % (seq, stderr))

    f = open(dot_fname)
    dot_lines = f.readlines()        
    f.close()

    mfe_struct = dot_lines[-1].strip('\n')

    os.remove(ct_fname)
    os.remove(dot_fname)

    return mfe_struct

def mfe_contrafold_(seq, T=37, version='2', constraint=None, param_file=None,DIRLOC=None,
    viterbi=False, probing_signal=None, probing_kws=None):
    """get MFE structure for Contrafold

    Args:
        seq (str): nucleic acid sequence
        T (float): temperature
        constraint (str): structure constraints
        motif (str): argument to vienna motif  
    Returns
        secondary structure dot-bracket string for MFE
    """
    if not version: version='2'

    if probing_signal is not None:
        fname = write_reactivity_file_contrafold(probing_signal, seq)
    else:
        fname = '%s.in' % filename()

    if DIRLOC is not None:
        LOC=DIRLOC
    elif version.startswith('2'):
        LOC=package_locs['contrafold_2']
    elif version.startswith('1'):
        LOC=package_locs['contrafold_1']
    else:
        raise RuntimeError('Error, Contrafold version %s not present' % version)

    command = ['%s/contrafold' % LOC, 'predict', fname]

    if probing_signal is not None:
        command = command + ['--evidence', '--params', package_locs['eternafold']+'/../parameters/EternaFoldParams_PLUS_POTENTIALS.v1', '--numdatasources','1', ]
        if probing_kws is not None:
            if 'kappa' in probing_kws.keys():
                command = command + ['--kappa', str(probing_kws['kappa']) ]
    else:
        if param_file is not None:
            command = command + ['--params', param_file]

    if viterbi:
        command.append('--viterbi')

    if constraint is not None:
        convert_dbn_to_contrafold_input(seq, constraint, fname)
        command.append('--constraints')
    else:
        if probing_signal is None:
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

def mfe_linearfold_(seq, bpps=False, package='contrafold', beam_size=100, return_dG_MFE=False):
    
    seqfile = write([seq])

    LOC = package_locs['linearfold']

    if bpps:

        pf_only = 0
    else:
        pf_only = 1

    # args:  beamsize, is_sharpturn, is_verbose, is_eval, is_constraints]
    #Todo: implement constraint input
    command=['echo %s | %s/linearfold_%s' % (seq, LOC, package[0]), str(beam_size), '0', '0', '0']
    if DEBUG: print(' '.join(command))
    p = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)

    stdout, stderr = p.communicate()

    if DEBUG:
        print('stdout')
        print(stdout)
        print('stderr')
        print(stderr)

    if p.returncode:
        raise Exception('LinearFold failed: on %s\n%s' % (seq, stderr))


    # linearfold returns two different things depending on which package
    struct = stdout.decode('utf-8').split('\n')[1].split(' ')[0]

    os.remove(seqfile)

    if return_dG_MFE:

        dG_mfe = float(stdout.decode('utf-8').split('\n')[1].split(' ')[1][1:-1])

        if package.lower() != 'vienna':
            dG_mfe *= -1

        return struct, dG_mfe

    else:
        return struct


