from arnie.mea.mea import *
import numpy as np
from glob import glob
import argparse
import sys, os

def predict_MEA_structures(matrix_list, gamma_min=-7, gamma_max=7, verbose=False, metric='mcc', output_dir='MEA_output'):
    '''Estimate maximum expected pseudoaccuracy structures per Hamada et al. BMC Bioinf 2010 11:586.
    
    Note: Files in matrix_dir and true_structs need to have the same names corresponding to their same constructs, but suffixes don't matter.
        
    Inputs:

    matrix_dir: list of NxN base pair probability matrices.
    gamma_min, gamma_max: min/max log_2(gamma) value used, defaults are -7 and 7.
    metric: keyword-based, which metric to use to select structure. Options are 'sen', 'ppv', 'mcc', 'fscore'.
    verbose: print output or not (for command line use)

    Outputs:
    List of predicted structures (in dbn format) at each gamma.

    '''

    metric_ind = ['sen', 'ppv', 'mcc', 'fscore'].index(metric)

    if len(matrix_list) == 0:
        raise ValueError('No matrix files found!')

    matrices = [np.loadtxt(x) for x in matrix_list]
    pdb_indices = [os.path.basename(x).split('.')[0] for x in matrix_list]

    n_constructs = len(matrices)

    gamma_vals = [x for x in range(gamma_min, gamma_max)]
    best_metric_values, best_gammas, best_structs,best_metrics = [],[],[],[]

    metrics_across_gammas = {k:[] for k in gamma_vals}

    if verbose: print('\nmetric\tpdb_ind\tbest_log2g\tbest_metric_value\tbest_struct')

    for i, matrix in enumerate(matrices):

        running_best_metrics = []
        running_best_value = 0
        running_best_gamma = -101
        running_best_struct = ''

        for g in gamma_vals:

            mea_cls = MEA(matrix, gamma=2**g)

            metrics = mea_cls.score_expected() #sen, ppv, mcc, fscore
            metrics_across_gammas[g].append(metrics)

            if metrics[metric_ind] > running_best_value:
                running_best_value = metrics[metric_ind]
                running_best_metrics = metrics
                running_best_gamma = g
                running_best_struct = mea_cls.structure

        best_metrics.append(running_best_metrics)
        best_metric_values.append(running_best_value)
        best_gammas.append(running_best_gamma)
        best_structs.append(running_best_struct)

        if verbose: print("%s\t%s\t%d\t%.3f\t%s" % (metric, pdb_indices[i], running_best_gamma, running_best_value, running_best_struct))

    # print('Avg metrics across gamma vals')

    print('\t\tlog2(g)\tsen\tppv\tmcc\tfscore')

    for g in gamma_vals:

        [sen, ppv, mcc, fscore] = np.mean(metrics_across_gammas[g], axis=0)
        print('gamma_avg\t%d\t%.3f\t%.3f\t%.3f\t%.3f' % (g, sen, ppv, mcc, fscore))

    # print('Best avg metrics using individual gammas')
    [sen, ppv, mcc, fscore] = np.mean(np.array(best_metrics), axis=0)

    print('gamma_best\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f' % (np.mean(best_gammas), sen, ppv, mcc, fscore))

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for struct, ind in list(zip(best_structs, pdb_indices)):
        if os.path.exists('%s/%s.dbn' % (output_dir, ind)):
            print('NB: overwriting existing predicted structure')
        with open('%s/%s.dbn' % (output_dir, ind), 'w') as f:
            f.write(struct)

    return best_structs
    
def score_against_true_structs(pred_struct_list, true_struct_list, verbose=False, weight_by_n_bps=False):
    '''Score maximum expected pseudoaccuracy structures against provided 3D structures.  
    
    Note: Files in matrix_dir and true_structs need to have the same names corresponding
        to their same constructs, but suffixes don't matter.
        
    Inputs:

    pred_struct_list: list of predicted structures.
    true_structs: list of NxN true structure base pair matrices. Can be
        symmetric matrices or not; upper triangle is taken.
    verbose: print output or not (for command line use)

    Outputs:

    SEN: TP/(TP+FN), library keyed by gamma values used.
    PPV: TP/(TP+FP), "
    MCC: Mathews correlation coefficient
    Fscore: 2*TP/(2*TP + FP + FN)    

    '''
    pred_structs, true_structs = [], []

    if len(pred_struct_list) == 0:
        raise ValueError('No predicted structure files found!')

    if len(true_struct_list) == 0:
        raise ValueError('No ground truth structure files found!')

    for x in pred_struct_list:
        for s in true_struct_list:
            if os.path.basename(x).split('.')[0] in s:

                pstruct = load_matrix_or_dbn(x)
                pred_structs.append(pstruct)

                struct = load_matrix_or_dbn(s)
                true_structs.append(struct)

    assert len(pred_structs) == len(true_structs)

    tally, ptl_sen, ptl_ppv, ptl_mcc, ptl_fscore = 0, 0, 0, 0, 0

    pdb_indices = [os.path.basename(x).split('.')[0] for x in pred_struct_list]
    
    for i in range(len(pred_structs)):
        
        sen, ppv, mcc, fscore, N = score_ground_truth(pred_structs[i], true_structs[i])
        print('Score:\t%s\t%.3f\t%.3f\t%.3f\t%.3f' % (pdb_indices[i], sen, ppv, mcc, fscore))

        if weight_by_n_bps:
            ptl_sen += sen*N
            ptl_ppv += ppv*N
            ptl_mcc += mcc*N
            ptl_fscore += fscore*N
            tally += N

        else:
            ptl_sen += sen
            ptl_ppv += ppv
            ptl_mcc += mcc
            ptl_fscore += fscore
            tally += 1

    mean_sen = ptl_sen/tally
    mean_ppv = ptl_ppv/tally
    mean_mcc = ptl_mcc/tally
    mean_fscore = ptl_fscore/tally

    print("Avg:\tsen\tppv\tmcc\tfscore\n\t%.3f\t%.3f\t%.3f\t%.3f" % (mean_sen, mean_ppv, mean_mcc, mean_fscore))

    return mean_sen, mean_ppv, mean_mcc, mean_fscore

if __name__ == '__main__':

    parser=argparse.ArgumentParser(
        description='''Estimate maximum expected pseudoaccuracy structures per Hamada et al. BMC Bioinf 2010 11:586 and\
        score against a ground truth dataset.\n

         Input format: Base pair probability matrices (specified in --bp_matrices) need to have same base names
          as structures (specified in --true_structs, and can be either dbn strings or NxN matrices),
           but the extensions for both types don't matter.''')

    parser.add_argument('--bp_matrices','-p', nargs='+', 
        help='path to NxN matrices of bp probabilities, i.e. `contrafold/*.bpps`.')

    parser.add_argument('--output_dir', '-o', 
        help="Path to output of predicted MEA structures. Default is `MEA_output`.", default = 'MEA_output')

    parser.add_argument('--true_structs','-s', nargs='+', 
        help='Optional: path to true structures, i.e. `rnaview/*.struct`. These can be dbn structures or NxN matrices.', default=None)

    parser.add_argument('--metric', default='mcc',
        help='Accuracy metric, options are `mcc`, `fscore`, `ppv`, or `sen`. Default is `mcc`.')

    parser.add_argument('--gamma_min',type=int, default=-7, help='Min value for log_2(gamma), default is -7')
    parser.add_argument('--gamma_max',type=int, default=7, help='Max value for log_2(gamma), default is 7')

    parser.add_argument('--weight_by_n_bps', dest='weight_by_n_bps', action='store_true', 
        help='For scoring to true structures, weight accuracy over dataset by number of bps.\
         If flag not included, equal weight across constructs.')

    parser.add_argument('--verbose', dest='verbose', action='store_true')
    parser.add_argument('--score_truth_only', dest='score_truth_only', action='store_true',
        help='Use if MEA structures already generated and only scoring to ground truth dataset.')

    #print help and exit if no args
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
        
    args = parser.parse_args()

    #if args.true_structs:
        #assert len(args.bp_matrices) == len(args.true_structs)

    if args.verbose:
        print('\nRNA MEA STRUCTURE PREDICTION')
        print('Number of structures: %d' % len(args.bp_matrices))
        print('Path to first base pair matrix: %s' % args.bp_matrices[0])
        if args.true_structs:
            print('Path to first true struct: %s' % args.true_structs[0])
        print('\nScanning gamma for MEA structure prediction:')

    if not args.score_truth_only:
        predict_MEA_structures(args.bp_matrices, gamma_min = args.gamma_min, gamma_max = args.gamma_max, verbose=args.verbose, metric = args.metric, output_dir = args.output_dir)

    if args.true_structs:
        if args.verbose: print('\nScoring provided true structures against maximum expected pseudoaccuracy structures:')
        score_against_true_structs(glob('%s/*' % args.output_dir), args.true_structs, verbose=args.verbose, weight_by_n_bps=args.weight_by_n_bps)
