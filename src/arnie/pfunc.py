import os, re, sys, shutil
import subprocess as sp
import random, string
import numpy as np
from .utils import *


# load package locations from yaml file, watch! global dict
package_locs = load_package_locations()

def pfunc(seq, package='vienna_2', T=37,
    constraint=None, motif=None, linear=False,
    dangles=True, noncanonical=False, pseudo=False, dna=False, DIRLOC=None,
    bpps=False, param_file=None, coaxial=True, reweight=None,
    return_free_energy = False, beam_size=100, DEBUG=False, threshknot=False,
    probing_signal=None, probing_kws = None):
    ''' Compute partition function for RNA sequence.

        Args:
        seq (str): nucleic acid sequence
        T (float): temperature (Celsius)
        constraint (str): structure constraints
        motif (str): argument to vienna motif
        linear (bool): call LinearPartition to estimate Z in Vienna or Contrafold
        pseudo (bool): nupack only, make prediction with pseudoknots
        dna (bool): nupack only, make prediction for DNA
        dangles (bool): dangles or not, specifiable for vienna, nupack
        coaxial (bool): coaxial stacking or not, specifiable for rnastructure, vfold
        noncanonical(bool): include noncanonical pairs or not (for contrafold, RNAstructure (Cyclefold))
        beam_size (int): beam size option for LinearPartition.
        threshknot (bool): call threshknot to predict pseudoknots (for contrafold, using LinearPartition)

        Possible packages:
        'vienna_2', 'vienna_1','contrafold_1','contrafold_2','nupack_95','nupack_99','rnasoft_2007','rnasoft_1999','rnastructure','vfold_0','vfold_1'

    Returns
        float: free energy
    '''

    if DEBUG: load_package_locations(DEBUG=True)

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
            print('Warning: LinearPartition only implemented for vienna, contrafold, eternafold.')

    # if pkg=='eternafold' and not 'eternafoldparams' in package_locs.keys():
    #     raise RuntimeError('Error: need to set path to EternaFold params to use eternafold hotkey.')


    if pseudo and pkg !='nupack':
        raise ValueError('pseudo only for use with nupack')

    if threshknot:
        if pkg !='contrafold':
            print('Warning: ThreshKnot only implemented for Contrafold')
        elif not linear:
            print('Warning: ThreshKnot only implemented with LinearPartition')

    if pkg=='vienna':
        if linear:
            Z, tmp_file = pfunc_linearpartition_(seq, package='vienna',bpps=bpps, beam_size=beam_size,
                return_free_energy=return_free_energy, DEBUG=DEBUG)

        else:
            Z, tmp_file = pfunc_vienna_(seq, version=version, T=T, dangles=dangles,
             constraint=constraint, motif=motif, bpps=bpps, param_file=param_file,
             reweight=reweight, return_free_energy=return_free_energy, DEBUG=DEBUG, probing_signal=probing_signal, probing_kws = probing_kws,)

    elif pkg=='contrafold':
        if linear:
            Z, tmp_file = pfunc_linearpartition_(seq, package='contrafold', bpps=bpps, beam_size=beam_size,
                return_free_energy=return_free_energy, DEBUG=DEBUG, threshknot=threshknot)
        else:
            Z, tmp_file = pfunc_contrafold_(seq, version=version, T=T,
                constraint=constraint, bpps=bpps, param_file=param_file, DIRLOC=DIRLOC,
                return_free_energy=return_free_energy, DEBUG=DEBUG)

    elif pkg=='rnastructure':
        Z, tmp_file = pfunc_rnastructure_(seq, version=version, T=T, coaxial=coaxial,
            constraint=constraint, bpps=bpps, return_free_energy=return_free_energy, DEBUG=DEBUG)

    elif pkg=='rnasoft':
        if constraint is not None:
            print("ERROR: RNAsoft is unable to handle constraints for calculating \
                partition functions, returning unconstrained Z.")

        Z, tmp_file = pfunc_rnasoft_(seq, version=version, T=T, constraint=constraint,
         bpps=bpps,return_free_energy=return_free_energy, DEBUG=DEBUG)

    elif pkg=='nupack':
        Z, tmp_file = pfunc_nupack_(seq, version=version, dangles=dangles, T=T, pseudo=pseudo, dna=dna, constraint=constraint,

            return_free_energy=return_free_energy, DEBUG=DEBUG)

    elif pkg=='vfold':
        Z, tmp_file = pfunc_vfold_(seq, version=version, T=T, coaxial=coaxial, DEBUG=DEBUG)

    elif pkg=='eternafold':
        if linear:
            Z, tmp_file = pfunc_linearpartition_(seq, package='eternafold', bpps=bpps, beam_size=beam_size,
                return_free_energy=return_free_energy, DEBUG=DEBUG)
        else:
            # Using contrafold code and eternafold params
            if 'eternafoldparams' in package_locs.keys() and 'eternafold' not in package_locs.keys():
                Z, tmp_file = pfunc_contrafold_(seq, T=T, constraint=constraint, 
                    bpps=bpps, param_file=package_locs['eternafoldparams'], DIRLOC=DIRLOC, return_free_energy=return_free_energy, DEBUG=DEBUG)

            elif 'eternafold' in package_locs.keys() and param_file is None:
                #Using eternafold code and params in eternafold codebase
                efold_param_file = os.environ['ETERNAFOLD_PARAMETERS'] if os.environ.get('ETERNAFOLD_PARAMETERS') else package_locs['eternafold']+'/../parameters/EternaFoldParams.v1'
                if not os.path.exists(efold_param_file):
                   raise RuntimeError('Error: Parameters not found at %s' % efold_param_file)
                else:
                    Z, tmp_file = pfunc_contrafold_(seq, T=T, constraint=constraint, 
                        bpps=bpps, param_file=efold_param_file, DIRLOC= package_locs['eternafold'], return_free_energy=return_free_energy, DEBUG=DEBUG)
            elif 'eternafold' in package_locs.keys() and param_file is not None:
                    Z, tmp_file = pfunc_contrafold_(seq, T=T, constraint=constraint, 
                        bpps=bpps, param_file=param_file, DIRLOC= package_locs['eternafold'], 
                        probing_kws=probing_kws, probing_signal = probing_signal, return_free_energy=return_free_energy, DEBUG=DEBUG)                

    else:
        raise ValueError('package %s not understood.' % package)

    if bpps:
        return Z, tmp_file

    else:
        if tmp_file:
            try:
                os.remove(tmp_file)
            except:
                pass
        return Z

def pfunc_vienna_(seq, T=37, version='2', constraint=None, motif=None, param_file=None,
                dangles=True, bpps=False, reweight=None, return_free_energy=False, DEBUG=False, probing_signal=None, shapeMethod='W', probing_kws=None):
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

    command = ['%s/RNAfold' % LOC, '-p', '-T', str(T)]


    if probing_signal is not None:
        if probing_kws is None:
            probing_kws={}
        probing_file = run_RNAPVmin(probing_signal, seq, LOC, DEBUG, **probing_kws)
        command.append('--shape=%s' % probing_file)
        command.append('--shapeMethod=%s' % shapeMethod)


    if version.startswith('2'):

        command.append('--bppmThreshold=0.0000000001')
        output_id = local_rand_filename()
        output_dot_ps_file = "%s_0001_dp.ps" % output_id
        command.append('--id-prefix=%s' % output_id)
    else:
        output_dot_ps_file = 'dot.ps'

    if motif is not None:
        command.append("--motif=%s" % motif)

    if constraint is not None:
        fname = write([seq, constraint])
        command.append('-C')
        if version=='2':
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

    if stderr.decode('utf-8').startswith('WARNING: '):
        print(stderr)

    if p.returncode:
        raise Exception('RNAfold failed: on %s\n%s' % (seq, stderr))
    os.remove(fname)

    if version.startswith('2'):
        os.remove("%s_0001_ss.ps" % output_id)

    if 'omitting constraint' in stderr.decode('utf-8'):
        free_energy = np.inf # Impossible structure
    else:
        m = re.search('([,|\(\.\)\]\[\{\}]+)\s+\[\s*(-*[0-9]+\.[0-9]+)', stdout.decode('utf-8'))

        free_energy = float(m.group(2))
        if DEBUG: print('free_energy: ', free_energy)

    if return_free_energy:
        return free_energy, output_dot_ps_file
    else: # return Z
        return np.exp(-1*free_energy/(.0019899*(273+T))), output_dot_ps_file

def pfunc_contrafold_(seq, T=37, version='2', constraint=None, bpps=False,
         param_file=None, return_free_energy=False, DIRLOC=None, DEBUG=False, probing_signal=None, probing_kws=None):

    """get partition function structure representation and free energy

    Args:
        seq (str): nucleic acid sequence
        T (float): temperature
        constraint (str): structure constraints
        motif (str): argument to vienna motif

        DIRLOC: sets location of contrafold specifically (Useful if there's several EternaFold builds to compare.)
    Returns
        float: partition function
        Note: If the constraint is impossible then Z wil be equal to the Z unconstrained
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

    if not os.path.isdir(LOC):
        command = ['%s' % LOC, 'predict', fname]
    else:
        command = ['%s/contrafold' % LOC, 'predict', fname]

    if probing_signal is not None:
        command = command + ['--evidence', '--params', package_locs['eternafoldparams_PLUS_POTENTIALS'], '--numdatasources','1', ]
        if probing_kws is not None:
            if 'kappa' in probing_kws.keys():
                command = command + ['--kappa', str(probing_kws['kappa']) ]
    else:
        if param_file is not None:
            command = command + ['--params', param_file]

    if bpps:
        posterior_fname = '%s.posteriors' % filename()
        command = command + ['--posteriors', '0.0000000001', posterior_fname]

    else:
        command.append('--partition')


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

    if not bpps:
        logZ = float(stdout.decode('utf-8').rstrip().split()[-1])

        if return_free_energy:
            return -1*logZ, None
        else:
            return np.exp(logZ), None
    else:
        return 0, posterior_fname

def pfunc_rnasoft_(seq, version='99', T=37, constraint=None, bpps=False, return_free_energy=False, DEBUG=False):
    DIR = package_locs['rnasoft']

    if not version: version='blstar'

    #note for mfe will use simfold instead of simfold pf

    # supported versions: 07, 99, 99-no-dangles, BL-no-dangles, BLstar, LAM-CG, NOM-CG
    param_locs = {'07': '%s/params/CG_best_parameters_ISMB2007.txt' % DIR,
    '99': '%s/params/turner_parameters_fm363_constrdangles.txt' % DIR,
    '99-no-dangles': '%s/params/turner_parameters_fm363_dangles0.txt' % DIR,
    'bl-no-dangles': '%s/params/BL-no-dangles.txt' % DIR,
    'blstar': '%s/params/BLstar.txt' % DIR,
    'lam-cg': '%s/params/LAM-CG.txt' % DIR,
    'nom-cg': '%s/params/NOM-CG.txt' % DIR}

    command = ['%s/simfold_pf' % DIR, '-s', seq, '-p', param_locs[version]]

    if DEBUG: print(' '.join(command))
    p = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE)

    stdout, stderr = p.communicate()

    if DEBUG:
        print('stdout')
        print(stdout)
        print('stderr')
        print(stderr)

    bpps_fname = '%s.bpps' % filename()
    if bpps:
        with open(bpps_fname,'w') as f:
            for line in stdout.decode('utf-8').split('\n')[5:]:
                if not 'Glog' in line and len(line) > 1:
                    f.write(line+'\n')

    if p.returncode:
        raise Exception('RNAsoft partition failed: on %s\n%s' % (seq, stderr))

    Z = float(stdout.decode('utf-8').split('\n')[1].split()[-1])

    if return_free_energy:
        return -1*np.log(Z), bpps_fname
    else:
        return Z, bpps_fname

def pfunc_nupack_(seq, version='95', T=37, dangles=True, constraint=None, return_free_energy=False, pseudo=False, dna=False, DEBUG=False):

    if not version: version='95'
    nupack_materials={'95': 'rna1995', '99': 'rna1999', 'dna':'dna1998'}

    if dna: version='dna'

    DIR = package_locs['nupack']

    if dangles:
        dangle_option='some'
    else:
        dangle_option='none'

    if constraint is not None:
        if '.' in constraint:
            print('Warning: NUPACK does not do constrained folding . and x')
        seqfile = write([seq, constraint])
        command=['%s/energy' % DIR, '%s' % seqfile.replace('.in',''),'-T', str(T),
             '-material', nupack_materials[version], '-dangles', dangle_option]

    else:
        seqfile = write([seq])
        command=['%s/pfunc' % DIR, '%s' % seqfile.replace('.in',''),'-T', str(T),
             '-material', nupack_materials[version], '-dangles', dangle_option]

    if pseudo:
        command.append('--pseudo')
    if DEBUG: print(' '.join(command))
    p = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE)

    stdout, stderr = p.communicate()

    if DEBUG:
        print('stdout')
        print(stdout)
        print('stderr')
        print(stderr)

    if p.returncode:
        raise Exception('Nupack pfunc failed: on %s\n%s' % (seq, stderr))

    if constraint is not None:
        free_energy = float(stdout.decode('utf-8').split('\n')[-2])
        Z=np.exp(-1*free_energy/(.0019899*(273+T)))

    else:
        free_energy = float(stdout.decode('utf-8').split('\n')[-3])
        Z=float(stdout.decode('utf-8').split('\n')[-2])

    os.remove(seqfile)

    if return_free_energy:
        return free_energy, None
    else:
        return Z, None

def pfunc_rnastructure_(seq, version=None, T=37, constraint=None, coaxial=True,
                            bpps=False, return_free_energy=False, DEBUG=False):
    """get partition function structure representation and free energy

    Args:
        seq (str): nucleic acid sequence
        T (float): temperature
        constraint (str): structure constraints
        motif (str): argument to vienna motif
        coaxial (bool): Coaxial stacking or not (default True)
    Returns
        float: partition function
    """

    seqfile = write([seq])
    pfsfile = '%s.pfs' % filename()
    DIR = package_locs['rnastructure']
    command = ['%s/partition' % DIR, seqfile, pfsfile, '-T', str(T+273)]

    if not coaxial:
        command.extend(['--disablecoax'])

    if constraint is not None:
        fname = '%s.CON' % filename()
        #print(fname)
        convert_dbn_to_RNAstructure_input(seq, constraint, fname)
        command.extend(['--constraint', fname])

    if DEBUG: print(' '.join(command))
    p = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE)

    stdout, stderr = p.communicate()

    if DEBUG:
        print('stdout')
        print(stdout)
        print('stderr')
        print(stderr)

    os.remove(seqfile)
    if constraint is not None:
        os.remove(fname)
    
    if p.returncode:
        raise Exception('RNAstructure partition failed: on %s\n%s' % (seq, stderr))

    if not bpps:
        command = ['%s/EnsembleEnergy' % DIR, pfsfile]

        if DEBUG: print(' '.join(command))
        p = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE)

        stdout, stderr = p.communicate()

        if DEBUG:
            print('stdout')
            print(stdout)
            print('stderr')
            print(stderr)

        if p.returncode:
            raise Exception('RNAstructure EnsembleEnergy failed: on %s\n%s' % (seq, stderr))

        if DEBUG: print(stdout.decode('utf-8').split('\n')[3])
        free_energy = float(stdout.decode('utf-8').split('\n')[3].split(' ')[-2])

        if return_free_energy:
            return free_energy, pfsfile
        else:
            return np.exp(-1*free_energy/(.0019*(273+T))), pfsfile
    else:
        return 0, pfsfile

def pfunc_vfold_(seq, version='0', T=37, coaxial=True, bpps=False, DEBUG=False):
    #available versions: 0 for Turner 04 params, 1 for Mfold 2.3 params
    #for bpps
        # command = ['%s/Vfold2d_npk_mac.o %d %d %s %s %d' % (DIR, int(coaxial),\
        #  T, infile, outfile, int(version))]

    DIR = package_locs["vfold"]

    cwd = os.getcwd()
    os.chdir(DIR) #vfold precompiled binaries don't work being called from elsewhere

    if DEBUG: print(os.getcwd())

    seqfile = write([seq])

    if sys.platform=="linux":
        platform='linux'
    elif sys.platform=="darwin":
        platform='mac'
    elif sys.platform=="win32":
        platform='win'
    else:
        raise RuntimeError('Vfold has binaries for linux, macOS, and win')


    command = ['./VfoldThermal_npk_%s.o %d %d %d %s tmp %d; cat tmp; rm tmp' % (platform, int(coaxial), T, T, seqfile, int(version))]

    if DEBUG: print(' '.join(command))

    p = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)

    stdout, stderr = p.communicate()
    os.chdir(cwd)

    if DEBUG:
        print('stdout')
        print(stdout)
        print('stderr')
        print(stderr)
    if p.returncode:
        raise Exception('VfoldThermal_npk failed: on %s\n%s' % (seq, stderr))

    Z=float(stdout.decode('utf-8').split('\n')[-2].split()[1])

    os.remove(seqfile)
    return Z, None
    #output: take second field of last line for Z


def pfunc_linearpartition_(seq, bpps=False, package='contrafold', beam_size=100, return_free_energy=False, DEBUG=False, threshknot=False):
    LOC = package_locs['linearpartition']
    tmp_file = filename()
    tmp_command = filename()

    if bpps:
        pf_only = 0
    else:
        pf_only = 1

    if threshknot:
        TK = '-T'
    else:
        TK = '_'

    # args: beamsize, is_sharpturn, is_verbose, bpp_file, bpp_prefix, pf_only, bpp_cutoff,
    #forest_file, mea, gamma, TK, threshold, ThreshKnot_prefix, MEA_prefix, MEA_bpseq

    # args: beamsize, is_sharpturn, is_verbose, bpp_file, bpp_prefix, pf_only, bpp_cutoff,
    #forest_file, mea, gamma, TK, threshold, ThreshKnot_prefix, MEA_prefix, MEA_bpseq, shape_file_path

    #threshknot threshold set to default 0.3

    command=['echo %s | %s/linearpartition_%s' % (seq, LOC, package[0].lower()), str(beam_size),
     '0', '0', tmp_file, "''", str(pf_only), '0.000001', "''", "''", "''",'%s' % (TK), "''", "''", "''", "''", "''"]

    with open('%s.sh' % tmp_command,'w') as f:
        f.write(' '.join(command))

    if DEBUG: print(' '.join(command))

    meta_command = ['chmod +x %s.sh; %s.sh' % (tmp_command, tmp_command)]
    p = sp.Popen(meta_command, stdout=sp.PIPE, stderr=sp.PIPE,shell=True)

    stdout, stderr = p.communicate(input=str.encode(seq))

    if DEBUG:
        print('stdout')
        print(stdout)
        print('stderr')
        print(stderr)

    if p.returncode:
        raise Exception('LinearPartition failed: on %s\n%s' % (seq, stderr))

    os.remove("%s.sh" % tmp_command)
    # Note: the linearfold exec says this is free energy in kcal/mol.
    # this is still just cfold log Z

    # linearfold returns two different things depending on which package

    if bpps:
        return 0, tmp_file
    else:

        if package in ['contrafold','eternafold']:
            logZ=float(stderr.decode('utf-8').split(' ')[-1])
            #os.remove('_') # tmp file written by linearpartition forest
            if return_free_energy:
                return -1*logZ, None
            else:
                return np.exp(logZ), None

        elif package=='vienna':
            free_energy = float(stderr.decode('utf-8').split(' ')[-2])
            T=37
            #os.remove('_') # tmp file written by linearpartition forest
            if return_free_energy:
                return free_energy, None
            else:
                return np.exp(-1*free_energy/(.0019899*(273+T))), None
