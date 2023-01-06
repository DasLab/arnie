from Bio import SeqIO
import arnie
from arnie.pk_predictors import pk_predict
from arnie.pk_predictors import pk_predict_from_bpp
from arnie.utils import *
from arnie.utils import _group_into_non_conflicting_bp
import pandas as pd
import numpy as np
from arnie.mfe_bootstrap import mfe_bootstrap
from arnie.bpps import bpps
from arnie.mea.mea import MEA
from arnie.mea.mea_utils import *
import sys, math
import argparse
import os
import time

def get_seq(seq_filename):
    '''takes fasta file with viral genome, returns uppercase string'''
    record = SeqIO.read(seq_filename, "fasta")
    dna_seq = str(record.seq).upper()
    return dna_seq.replace("T", "U")

def normalize_shape(shape_reacs):
    '''takes a list of shape values as floats and np.nan objects and returns normalized shape data
    'box-plot' normalize:
    -- remove outliers, i.e., any values above 1.5 * interquartile range
    -- Find maximum value after filtering. If its smaller than the 95th percentile value, then take that instead.
    -- scalefactor is mean of top 10th percentile of values, but removing values above that filter.'''

    shape_reacs = np.array(shape_reacs)
    # Get rid of nan values for now
    nonan_shape_reacs = shape_reacs[~np.isnan(shape_reacs)]
    # Find max filter 1: 1.5 * Inter-Quartile Range # RCK added note to make clear we are only filtering high (above Q3) values
    sorted_shape = np.sort(nonan_shape_reacs)
    
    q1 = sorted_shape[int(0.25 * len(sorted_shape))]
    q3 = sorted_shape[int(0.75 * len(sorted_shape))]
    iq_range = abs(q3 - q1)
    # RCK this is strange to me, normally when you use the IQR method, outliers are 1.5 IQR above Q3 or 1.5 IQR below Q1, herre we are looking at 1.5 IQR above 0
    # if everythin within 1.5 * IQR do not need to filter, otherwise it will errori
    if sorted_shape[-1] > 1.5 * iq_range:
        filter1 = next(x for x, val in \
            enumerate(list(sorted_shape)) if val > 1.5 * iq_range)

        # Find max Filter 2: 95% value
        filter2 = int(0.95 * len(sorted_shape))

        # Get maximum filter value and fiter data
        filter_cutoff = sorted_shape[max(filter1, filter2)]
        sorted_shape = sorted_shape[sorted_shape < filter_cutoff]
    # Scalefactor: Mean of top 10th percentile of values
    top90 = sorted_shape[int(0.9 * len(sorted_shape))]
    scalefactor = np.mean(sorted_shape[sorted_shape >= top90]) # RCK made >= incase top 10th percentile is only the last value
    # Scale dataset
    return shape_reacs/scalefactor

def get_normalized_shape_data(shape_filename):
    '''takes text file of any shape data set and returns a normalized list of values in a list
    (containing np.nan objects)'''

    shape_reacs = np.loadtxt(shape_filename, delimiter ="\n")
    shape_reacs[shape_reacs == -999] = np.nan
    # normalize shape data
    normalized_shape_data = normalize_shape(shape_reacs).tolist()
    return normalized_shape_data

def get_sliding_windows(full_seq, step, window):
    '''takes seq of viral genome and returns a list of windows and their starting coordinates'''
    coords = list(range(0,len(full_seq)-window+1,step))
    windows = []
    for i in coords:
        new_window = full_seq[i:i+window]
        windows.append(new_window)
    return windows, coords

def get_structure(seq, coord, pk_predictor, window, bpp_package='contrafold', linear_partition=True):
    '''takes the list of windows, coords, and the desired predictor and outputs a dataframe with the predicted
    structures for every window and a column indicating whether or not it is a pseudoknot and whether or not
    it is "stable"'''
    if pk_predictor == 'threshknot':
        bpp = bpps(seq, package=bpp_package, linear=linear_partition)# RCK maybe make this and option to turn on DEBUG otherwise it just prints a bunch of stuff, DEBUG=True)
        dotbracket = pk_predict_from_bpp(bpp, heuristic='threshknot')
        # RCK added pk_predictor to tell difference of linear parition
        pk_predictor_str = f'{pk_predictor}_{bpp_package}'
        if linear_partition:
            pk_predictor_str += '_linearpartition'
    else:
        dotbracket = pk_predict(seq, pk_predictor)
        pk_predictor_str = pk_predictor
    return pk_predictor_str, int(coord), int(coord)+window, seq, dotbracket, is_PK(dotbracket)

def shapeknots_predict(seq, shape, track, coord, window):
    '''takes seq windows, windows of a single track of shape data, coords, and window length;
    returns two dataframes with predicted structures for the genome, one with shapeknots mfe structures
    and the other with probknot mea structures'''

    shapeknot_struct, bpp = mfe_bootstrap(seq=seq, num_bootstrap=100, shape_signal=shape, pk=True)
    probknot_struct = pk_predict_from_bpp(bpp, heuristic='threshknot')
    shape_struct_list = []
    shape_struct_list.append([track+'_shapeknots', int(coord), int(coord)+int(window), seq, shapeknot_struct, is_PK(shapeknot_struct)])
    shape_struct_list.append([track+'_shape_probknot', int(coord), int(coord)+int(window), seq, probknot_struct, is_PK(probknot_struct)])

    return shape_struct_list


def unstable_helices(dotbracket):
    '''takes dotbracket structure of a pseudoknot and returns True if the pseudoknotted helices have at least 2 bps'''
    # note - might update this function so that it saves the structure if at least some of the pseudoknotted helices
    ## have 2 base pairs, but doesn't automatically eliminate structures that have one helix with one bp
    bp_list = convert_dotbracket_to_bp_list(dotbracket, allow_pseudoknots=True)
    groups = _group_into_non_conflicting_bp(bp_list)
    for idx, pairs in enumerate(groups):
        if idx > 0:
            if len(pairs) == 1:
                return True
    # RCK let's chat about this function, it is not great practice for there
    # to be scenarios where a function just doesn't return anything
    # also I agree with you, I think it is fine to delete helices with len < 1
    # but this function isn't doing that.

    #GPN sounds great - this function may not even be necessary,
        #I think I can find another way to accomplish my real goal

def determine_outliers(data):
    '''takes a list of values and returns a list of True or False statements for each value
    stating whether or not they represent statistical outliers'''
    sorted_data = np.sort(data)
    q1 = sorted_data[int(0.25 * len(sorted_data))]
    q3 = sorted_data[int(0.75 * len(sorted_data))]
    iqr = q3-q1
    lower_bound = q1-(1.5*iqr)
    upper_bound = q3+(1.5*iqr)
    # RCK using numpt again you can just combine these statements
    outliers = (data < lower_bound) | (data > upper_bound)
    return outliers

def evaluate_L1_shape_score(s,shape):
    '''takes the structure in dotbracket format and a list of shape reactivity values of the same length
    and returns a score from 0-1 on the agreement of the structure and the shape data'''
    score = 0
    for c,react in zip(s,shape):
        if (c=="." and react>0.25) or (c!="." and react<0.5):
            score += 1
    return score/len(s)

def get_pk_bp_locs(dotbracket):
    '''takes a dotbracket structure and returns a list of bases in the string which are involved in
    pseudoknotted helices'''
    pk_bp_locs = []
    bp_list = convert_dotbracket_to_bp_list(dotbracket, allow_pseudoknots=True)
    groups = _group_into_non_conflicting_bp(bp_list)
    for helix in groups[1:]:

            for i in range(len(helix)):
                pk_bp_locs += helix[i]
    return sorted(pk_bp_locs)

#TODO: add a sed command to give each folder its own out and error folder and delete all after
#probably need to change the arguments to this function so that it inputs multiple sequences to the sbatch file
def get_struct_on_node(seq_windows_str, coords_str, template_sbatch, temp_folder, pk_predictors_str, window, bpp_package='contrafold',linear_partition=True):
    coords_list = coords_str.split(" ")
    sbatch_name = ''
    for i in range(len(coords_list)):
        sbatch_name += coords_list[i] + '_'
    os.system(f'cp {template_sbatch} {temp_folder}/{sbatch_name}.sbatch')
    file = open(f"{temp_folder}/{sbatch_name}.sbatch", "a")
    # RCK you had hard encoded you scratch folder, changed this to check the same path as current script 
    file.write(f"python {os.path.dirname(__file__)}/ViralKnots_single.py --seqs {seq_windows_str} --coords {coords_str} --temp_folder {temp_folder} --pk_predictors {pk_predictors_str} -w {window} --bpp_package {bpp_package} --linear_partition {linear_partition} \n")
    file.close()
    os.popen(f'sbatch {temp_folder}/{sbatch_name}.sbatch')

def get_shape_on_node(seq_windows_str, coords_str, shape_windows_str, template_sbatch, temp_folder, shape_data_sets_str, window):
    coords_list = coords_str.split(" ")
    sbatch_name = ''
    for i in range(len(coords_list)):
        # RCK added shape so does not overwrite pk_predict outputs
        sbatch_name += coords_list[i] + '_'
    os.system(f'cp {template_sbatch} {temp_folder}/{sbatch_name}_shape.sbatch')
    file = open(f"{temp_folder}/{sbatch_name}_shape.sbatch", "a")
    # RCK you had hard encoded you scratch folder, changed this to check the same path as current script 
    file.write(f"python {os.path.dirname(__file__)}/ViralKnots_single_shape.py --seqs {seq_windows_str} --coords {coords_str} --shape_windows {shape_windows_str} --temp_folder {temp_folder} --shape_data_sets {shape_data_sets_str} -w {window} \n")
    file.close()
    os.popen(f'sbatch {temp_folder}/{sbatch_name}_shape.sbatch')

def combine_struct_files(temp_folder):
    dfs = []
    for file in os.listdir(temp_folder):
        if file[-4:] == '.csv':
            df = pd.read_csv(f"{temp_folder}/{file}")
            dfs.append(df)
    return dfs

