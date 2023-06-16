import random
import numpy as np
from .mfe import mfe
from .utils import get_bpp_from_dbn
from .utils import filename
from .utils import load_package_locations
from os import remove

# load package locations from yaml file, watch! global dict
package_locs = load_package_locations()

def get_bootstrap_reac_file(reactivity):
    reac_file = '%s.SHAPE' % filename()
    range_arr = np.arange(1, len(reactivity) + 1)
    reac_arr = np.array(reactivity)
    shape_pos = np.array([range_arr, reac_arr]).T
    sample_idx = np.random.choice(len(reactivity), len(reactivity))
    shape_pos = shape_pos[sample_idx,:]

    f = open(reac_file, 'w')
    for cur_sample in shape_pos:
        pos, reactivity = cur_sample
        if reactivity > 0:
            f.write('%d %f\n' % (pos, reactivity))
    f.close()

    return reac_file

def mfe_bootstrap(seq, num_bootstrap, 
    package='rnastructure', T=37,
    constraint=None, shape_signal=None, dms_signal=None, pseudo=False):
    """
    Compute MFE structure (within package) for RNA sequence with bootstrapping on the SHAPE/DMS data.

        Args:
        seq (str): nucleic acid sequence
        T (float): temperature (Celsius)
        constraint (str): structure constraints
        shape_signal(list): list of normalized SHAPE reactivities, with negative values indicating no signal
        dms_signal(list): list of normalized DMS reactivities, with negative values indicating no signal
        pseudo: if True, will predict pseudoknots, but only with RNAstructure

        Possible packages: 
        'rnastructure'
        
    Returns
        string: MFE structure
        np array: Base-pair probability matrix from bootstrapping
    """
    if (shape_signal is None) and (dms_signal is None):
        raise ValueError("Bootstrapping only applies if you have reactivity data.")
    if package != 'rnastructure':
        raise ValueError("Bootstrapping only runs for now with RNAstructure")

    bpp_matrix = np.zeros((len(seq), len(seq)))

    mfe_struct = mfe(seq, package=package, T=T, constraint=constraint, 
        shape_signal=shape_signal, dms_signal=dms_signal, pseudo=pseudo)

    for bootstrap in range(num_bootstrap):
        shape_file = None
        dms_file = None
        
        if shape_signal is not None:
            shape_file = get_bootstrap_reac_file(shape_signal)
        if dms_signal is not None:
            dms_file = get_bootstrap_reac_file(dms_signal)

        cur_mfe_struct = mfe(seq, package=package, T=T, constraint=constraint, 
            shape_file=shape_file, dms_file=dms_file, pseudo=pseudo)
        bpp_matrix += get_bpp_from_dbn(cur_mfe_struct)

        if shape_signal is not None:
            remove(shape_file)
        if dms_signal is not None:
            remove(dms_file)

    return [mfe_struct, bpp_matrix/num_bootstrap]
