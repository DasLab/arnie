from ViralKnots_utils import *

def viral_knots(seq_filename, step, window, pk_predictors=[], pk_predict=False, bpp_package='contrafold', shape_data_folder=None, shape_data_sets=[], shapeknots=False,
                shape_rankings=False, spawn=False, template_sbatch=None, num_jobs=None, linear_partition=True):
    '''Takes viral genome and outputs csv with predicted structures, including pseudoknots
        Args:
            seq_filename - fasta file containing RNA sequence for viral genome
            step - desired number of nucleotides to slide each window
            window - desired size of window to divide genome into
            pk_predict - run pk predictions with algorithms other than shapeknots
            pk_predictors - list of desired predictors for use: threshknots, knotty, pknots, spotrna
            shapeknots - run pk predictions with shapeknots (must include shape data)
            shape_rankings - use shape reactivity values to score the pseudoknots (must include shape data)
            shape_data_folder - folder containing csv's with shape reactivity
            shape_data_sets - names of the reactivity files with a set of data the same length as the viral genome
                (can include nan values if full data for the genome is not available)
                (no .csv necessary)
                (note that this is designed to work with single or multiple tracks of shape data)'''
    ##first retrieve viral sequence and normalized shape data (if directed) and sort into windows
    seq = get_seq(seq_filename)
    seq_windows, coords = get_sliding_windows(seq, step=step, window=window)

    shape_sets = []
    #all_shape_windows is a list with all five shape sets subdivided into windows
    all_shape_windows = []
    if shapeknots or shape_rankings:
        for data_set in shape_data_sets:
            shape_data = get_normalized_shape_data(shape_data_folder+'/'+data_set+'.csv')
            # RCK prevent errors when shape data is not correct length
            assert len(seq) == len(shape_data), f'The length of sequence in the fasta file is {len(seq)}, but the length of shape data in {data_set} is {len(shape_data)}'
            shape_sets.append(shape_data)

        for shape_set in shape_sets:
            shape_windows, shape_coords = get_sliding_windows(shape_set, step=step, window=window)
            all_shape_windows.append(shape_windows)
    #run pk predictors, if directed
    dfs = []
    if spawn:
        temp_folder = get_random_folder(8)
        os.mkdir(temp_folder)
        # RCK I thought there was a bug in this and I couldn't read it so I rewrote it, turned out bug was elsewhere and of my making, but hey heres shorter code.
        num_tasks = 0
        if pk_predict:
            num_tasks += len(seq_windows)
        if shapeknots:
            num_tasks += len(seq_windows)+len(shape_data_sets)
        size_job = math.ceil(num_tasks/num_jobs)

    if pk_predict:
        if spawn:
            #now formatting to be list separated by spaces for input into ViralKnots_single
            pk_predictors_str = " ".join(pk_predictors)
            for i in range(0,len(seq_windows),size_job):
                seq_windows_str = " ".join([str(x) for x in seq_windows[i:i+size_job]])
                coords_str = " ".join([str(x) for x in coords[i:i+size_job]])
                get_struct_on_node(seq_windows_str, coords_str, template_sbatch, temp_folder, pk_predictors_str, window=window, bpp_package=bpp_package,linear_partition=linear_partition)

        else:
            for name in pk_predictors:
                struct_list = []
                for seq,coord in zip(seq_windows,coords):
                    struct_list.append(get_structure(seq, coord, pk_predictor=name, window=window, bpp_package=bpp_package,linear_partition=linear_partition))
                df = pd.DataFrame(struct_list,columns=["predictor","start","end","sequence", "struct", "pseudoknot"])
                dfs.append(df)
    #return coord, coord+window, seq, dotbracket, is_PK(dotbracket)
    #run shapeknots, if directed
    #TO DO: make shapeknots_predict logic looped through in this larger function for parallelization
    if shapeknots:
        if spawn:
            shape_data_sets_str = " ".join(shape_data_sets)
            for track, shape_set in enumerate(all_shape_windows):
                for i in range(0,len(seq_windows),size_job):
                    seq_windows_str = " ".join([str(x) for x in seq_windows[i:i+size_job]])
                    coords_str = " ".join([str(x) for x in coords[i:i+size_job]])
                    shape_windows_str = " ".join([",".join(map(str, x)) for x in shape_set[i:i+size_job]])
                    get_shape_on_node(seq_windows_str, coords_str, shape_windows_str, template_sbatch, temp_folder, shape_data_sets_str, window=window)
        else:
            for track, shape_set in enumerate(all_shape_windows):
                shapeknot_struct_list = []
                probknot_struct_list = []
                for seq,shape,coord in zip(seq_windows, shape_set, coords):
                    shape_df, prob_df = shapeknots_predict(seq, shape, shape_data_sets[track], coord,
                                                           window=window)
                    dfs.append(shape_df)
                    dfs.append(prob_df)

    if spawn:
        # RCK there won't always be this many jobs spawned so you want the min no?
        # while (2*(len(seq_windows)*len(pk_predictors)+2*(len(seq_windows)*len(shape_data_sets)))) > len(os.listdir(temp_folder)):
        # all results are in once there is a sbatch script and output csv written for each job.
        # The number of jobs run is the user specified, num_jobs, if there were more tasks, but otherwise num_tasks
        while 2*min(num_tasks,num_jobs) > len(os.listdir(temp_folder)):
            time.sleep(5) # check every 5 seconds
        dfs = combine_struct_files(temp_folder)
        files = os.listdir(temp_folder)
        for f in files:
            os.remove(f"{temp_folder}/{f}")
        os.rmdir(temp_folder)

    df1 = pd.concat(dfs)

    #first generate subsets of smaller dataframes with a single window
    #then calculate the sensitivity and ppv of all structures predicted for that window
    df2s = []
    for i,coord in enumerate(coords):

        df2 = df1.loc[df1['start'] == coord].copy() # RCK to get rid of warning, since its a copy now also got rid of the .loc fix below (x8)
        dotbrackets = df2['struct'].to_list()

        F1_scores_for_window = []
        F1_scores_for_pks_in_window = []
        
        if len(pk_predictors)==1:
            df2['F1_score'] = np.nan
            df2['F1_outlier'] = np.nan
            df2['F1_for_pk_bps'] = np.nan
            df2['F1_for_pk_bps_outlier'] = np.nan
        
        #TO DO: if is not a pseudoknot, should give np.nan value for F1_scores_for_pk
        else:
            for idx, dotbracket1 in enumerate(dotbrackets):
                F1_scores_for_struct = []
                F1_scores_for_pk_struct = []
                for idx2, dotbracket2 in enumerate(dotbrackets):
                    if idx != idx2:

                        F1_scores_for_struct.append(compare_structures_to_natives([dotbracket1], [dotbracket2],
                                                                     comparison='basepairs', metric='F1_score'))

                        F1_scores_for_pk_struct.append(compare_structures_to_natives([dotbracket1], [dotbracket2],
                                                                     comparison='PK_basepairs', metric='F1_score'))

                F1_scores_for_window.append((sum(F1_scores_for_struct))/(len(F1_scores_for_struct)))

                F1_scores_for_pks_in_window.append((sum(F1_scores_for_pk_struct))/(len(F1_scores_for_pk_struct)))

            #add all F1 scores and outlier determinations to the dataframe as new columns
            df2['F1_score'] = F1_scores_for_window
            df2['F1_outlier'] = determine_outliers(F1_scores_for_window)
            df2['F1_for_pk_bps'] = F1_scores_for_pks_in_window
            df2['F1_for_pk_bps_outlier'] = determine_outliers(F1_scores_for_pks_in_window)
        
        if shape_rankings:
            shape_sets_for_window = []
            for track in all_shape_windows:
                shape_sets_for_window.append(track[i])

            shape_scores_for_all_structs_in_window = []
            for struct in dotbrackets:
                struct_scores = []
                for shape_set in shape_sets_for_window:
                    struct_scores.append(evaluate_L1_shape_score(struct, shape_set))
                shape_scores_for_all_structs_in_window.append(sum(struct_scores)/len(struct_scores))
            df2.loc[:,'shape_score'] = shape_scores_for_all_structs_in_window

            #calculate the shape agreement of the pseudoknotted helices in every structure

            pk_bp_shape_scores_for_all_structs_in_window = []
            for struct in dotbrackets:
                pk_bp_locs = get_pk_bp_locs(struct)
                abbv_dotbracket = ''
                #all_abbv_shape_sets_for_struct is a set of five shape sets (one for each dataset)
                 #corresponding with only pk bps]
                for loc in pk_bp_locs:
                    abbv_dotbracket += struct[loc]
                all_abbv_shape_sets_for_struct = []
                for shape_set in shape_sets_for_window:
                    abbv_shape_set_for_struct = []
                    for loc in pk_bp_locs:
                        abbv_shape_set_for_struct.append(shape_set[loc])
                    all_abbv_shape_sets_for_struct.append(abbv_shape_set_for_struct)
                all_pk_bp_shape_scores_for_struct = []
                for abbv_shape in all_abbv_shape_sets_for_struct:
                    if len(abbv_shape) != 0:
                        all_pk_bp_shape_scores_for_struct.append(evaluate_L1_shape_score(abbv_dotbracket, abbv_shape))
                    else:
                        all_pk_bp_shape_scores_for_struct.append(0)

                if len(all_pk_bp_shape_scores_for_struct) != 0:
                    pk_bp_shape_scores_for_all_structs_in_window.append(sum(all_pk_bp_shape_scores_for_struct)/
                                                                len(all_pk_bp_shape_scores_for_struct))
                else:
                    pk_bp_shape_scores_for_all_structs_in_window.append(0)

            df2.loc[:,'pk_bp_shape_score'] = pk_bp_shape_scores_for_all_structs_in_window

        df2s.append(df2)

    df3 = pd.concat(df2s)
    df3.to_csv('output.csv', index=False)
    return None

if __name__=='__main__':

    parser=argparse.ArgumentParser()
    
    # RCK may be useful to add --help description for these inputs, example is below
    parser.add_argument('--seq_filename', '-s', type=str, required=True, help="Sequence file in fasta format") 
    parser.add_argument('--step', type=int, required=True)
    parser.add_argument('--window', '-w', type=int, required=True)
    parser.add_argument('--pk_predictors', default=[], nargs='+') 
    parser.add_argument('--pk_predict', action='store_true') # RCK changed so instead of saying true or false just needs to put '--pk_predict' if true and nothing if not
    parser.add_argument('--shapeknots', action='store_true') # RCK same as pk_predict
    parser.add_argument('--shape_data_folder', default='.', type=str) # RCK I am lazy so just set current directory to default
    parser.add_argument('--shape_data_sets', default=[], nargs='+', help='name of text files in shape_data_folder that has shape reactivity values for each nucleotide, 1 per line, and -999 if unkown, msut have a .csv extension but .csv should not be included in the name given') # RCK added help because I got tripped up on the format of this argument
    parser.add_argument('--shape_rankings', action='store_true') # RCK same as pk_predict
    parser.add_argument('--bpp_package', default='contrafold_2', type=str)
    parser.add_argument('--linear_partition', action='store_true') # RCK added
    parser.add_argument('--spawn', action='store_true') # RCK same as pk_predict
    parser.add_argument('--template_sbatch', default=None, type=str)
    #parser.add_argument('--temp_folder', default=None, type=str) # RCK there was already an arnie function for randomly making this tmp folder so I just used that?
    parser.add_argument('--num_jobs', default=None, type=int)

    args=parser.parse_args()
    
    # RCK added checks to get reasonable errors
    assert args.pk_predict or args.shapeknots, 'have to run --pk_predict and/or --shapeknots'
    assert not (args.pk_predict and args.pk_predictors==[]), 'if running pk_predict, need to specify --pk_predictors'
    assert not (args.spawn and ((args.template_sbatch is None) or (args.num_jobs is None))), 'if spawning jobs need to specify a --template_sbatch and --num_jobs'
    assert not (args.shapeknots and args.shape_data_sets==[]), 'if running shapeknots need to provide --shape_data_sets and --shape_data_folder if shape file are eleswhere'
    assert not (args.shape_rankings and args.shape_data_sets==[]), 'if running shaperankings need to provide --shape_data_sets and --shape_data_folder if shape file are eleswhere'
    assert not (args.spawn and args.shapeknots and args.pk_predict and args.num_jobs==1), 'to spawn with shapeknots and pk_predict need at least --num_jobs 2'
    
    # RCK added to check package compatability
    package_str = "Running PK-predictors:"
    for package in args.pk_predictors:
        if package == "threshknot":
            if args.linear_partition:
                assert args.bpp_package in ['vienna_2', 'contrafold_2', 'eternafold','contrafold','vienna'], "LinearPartition only implemented for vienna_2, contrafold_2, eternafold, vienna, contrafold."
                package_str += f' threshknot with base-pair-probaility matrix from {args.bpp_package} with linear partition,'
            else:
                package_str += f' threshknot with base-pair-probaility matrix from {args.bpp_package},'
        else:
            package_str += f' {package},'
    print(package_str[:-1])
    
    viral_knots(args.seq_filename, args.step, args.window, args.pk_predictors, args.pk_predict, args.bpp_package, args.shape_data_folder, args.shape_data_sets, args.shapeknots, args.shape_rankings, args.spawn, args.template_sbatch, args.num_jobs, args.linear_partition)

