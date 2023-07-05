import subprocess as sp
from arnie.utils import *
import glob
from os import getcwd, chdir, remove, mkdir, rmdir, path
from scipy.optimize import linear_sum_assignment


# TODO script all previous investigations
# TODO Debug modes to print output and err to help with install issues
# TODO pk_predict options +

package_locs = load_package_locations()


def pk_predict(seq, predictor,
               model="default", param="parameters_DP03.txt",
               refinement=1, t1="auto", t2='auto',
               cpu=32):
    '''

    ipknot options:
        model: one of ["LinearPartition-C","LinearPartition-V","Boltzmann","ViennaRNA","CONTRAfold","NUPACK"]
        t1: probability threshold level 1 
        t2: probability threshold level 2
        refinement: number of times for refinment

    hotknots options:
        model: one of ["CC","RE","DP"]
        param: one of ["parameters_CC06.txt","parameters_CC09.txt","parameters_DP03.txt","parameters_DP09.txt"]

    spotrna options:
        cpu: number cpu threads

    e2efold options:
        ???

    nupack options:
        ????

    '''
    if predictor not in ["hotknots", "ipknot", "knotty", "spotrna", "e2efold", "pknots","spotrna2","nupack"]:
        raise ValueError('Only hotknots,ipknot,knotty,spotrna,spotrna2,e2efold,pknots,nupack implemented.')
    if predictor == "spotrna":
        return _run_spotrna(seq, cpu=cpu)[0]
    elif predictor == "spotrna2":
        return _run_spotrna2(seq)[0]
    elif predictor == "e2efold":
        return _e2efold(seq)
    elif predictor == "pknots":
        return _pknots(seq)
    elif predictor == "knotty":
        return _knotty_mfe(seq)
    elif predictor == "hotknots":
        if model == "default":
            model = "DP"
        if model not in ["CC", "RE", "DP"]:
            raise ValueError('Only CC, RE, DP model implemented for hotknots.')
        if param not in ["parameters_CC06.txt", "parameters_CC09.txt", "parameters_DP03.txt", "parameters_DP09.txt"]:
            raise ValueError('Only parameters_CC06.txt, parameters_CC09.txt, parameters_DP03.txt, parameters_DP09.txt parameters implemented for hotknots.')
        return _run_hotknots(seq, model=model, param=param)[0][0]
    elif predictor == "ipknot":
        if model == "default":
            model = "LinearPartition-C"
        if model not in ["LinearPartition-C", "LinearPartition-V", "Boltzmann", "ViennaRNA", "CONTRAfold", "NUPACK"]:
            raise ValueError('Only LinearPartition-C, LinearPartition-V, Boltzmann, ViennaRNA, CONTRAfold, NUPACK model implemented for ipknot.')
        return _ipknot_mfe(seq, model=model, refinement=refinement, t1=t1, t2=t2)
    elif predictor == "nupack":
        return _nupack_mfe_pk(seq)


def pk_predict_from_bpp(bpp, heuristic="hungarian", theta=None, allowed_buldge_len=0, min_len_helix=2,
                        exp=1, sigmoid_slope_factor=None, prob_to_0_threshold_prior=0, prob_to_1_threshold_prior=1, ln=False, add_p_unpaired=True,
                        max_iter=1):
    '''
    threshknot options:
        theta
        max_iter
        allowed_buldge_len
        min_len_helix

    hungarian options:
        add_p_unpaired
        theta (aka prob_to_0_threshold_post)
        prob_to_0_threshold_prior
        prob_to_1_threshold_prior
        exp
        sigmoid_slope_factor
        ln
        allowed_buldge_len
        min_len_helix
    '''

    if heuristic not in ["threshknot", "hungarian"]:
        raise ValueError('Only threshknot and hunagrian heuristics implemented.')

    if heuristic == "threshknot":
        if theta is None:
            theta = 0.3
        return _threshknot(bpp, theta=theta, max_iter=max_iter, allowed_buldge_len=allowed_buldge_len, min_len_helix=min_len_helix)[0]
    elif heuristic == "hungarian":
        if theta is None:
            theta = 0.0
        return _hungarian(bpp, exp=1, sigmoid_slope_factor=sigmoid_slope_factor, prob_to_0_threshold_prior=prob_to_0_threshold_prior,
                          prob_to_1_threshold_prior=prob_to_1_threshold_prior, theta=theta, ln=ln, add_p_unpaired=add_p_unpaired,
                          allowed_buldge_len=allowed_buldge_len, min_len_helix=min_len_helix)[0]


def _hungarian(bpp, exp=1, sigmoid_slope_factor=None, prob_to_0_threshold_prior=0,
               prob_to_1_threshold_prior=1, theta=0, ln=False, add_p_unpaired=True,
               allowed_buldge_len=0, min_len_helix=2):

    bpp_orig = bpp.copy()

    if add_p_unpaired:
        p_unpaired = 1 - np.sum(bpp, axis=0)
        for i, punp in enumerate(p_unpaired):
            bpp[i, i] = punp

    # apply prob_to_0 threshold and prob_to_1 threshold
    bpp = np.where(bpp < prob_to_0_threshold_prior, 0, bpp)
    bpp = np.where(bpp > prob_to_1_threshold_prior, 1, bpp)

    # aply exponential. On second thought this is likely not as helpful as sigmoid since
    # * for 0 < exp < 1 lower probs will increase more than higher ones (seems undesirable)
    # * for exp > 1 all probs will decrease, which seems undesirable (but at least lower probs decrease more than higher ones)
    bpp = np.power(bpp, exp)

    # apply log which follows botlzamann where -ln(P) porportional to Energy
    if ln:
        bpp = np.log(bpp)

    bpp = np.where(np.isneginf(bpp), -1e10, bpp)
    bpp = np.where(np.isposinf(bpp), 1e10, bpp)

    # apply sigmoid modified by slope factor
    if sigmoid_slope_factor is not None and np.any(bpp):
        bpp = _sigmoid(bpp, slope_factor=sigmoid_slope_factor)

        # should think about order of above functions and possibly normalize again here

        # run hungarian algorithm to find base pairs
    _, row_pairs = linear_sum_assignment(-bpp)
    bp_list = []
    for col, row in enumerate(row_pairs):
        # if bpp_orig[col, row] != bpp[col, row]:
        #    print(col, row, bpp_orig[col, row], bpp[col, row])
        if bpp_orig[col, row] > theta and col < row:
            bp_list.append([col, row])

    structure = convert_bp_list_to_dotbracket(bp_list, bpp.shape[0])
    structure = post_process_struct(structure, allowed_buldge_len, min_len_helix)
    bp_list = convert_dotbracket_to_bp_list(structure, allow_pseudoknots=True)

    return structure, bp_list


def _sigmoid(x, slope_factor=0.5):
    # normalize to [-1, 1]
    numerator = (x - x.min()) * 2.0
    denominator = x.max() - x.min()
    #print(numerator, denominator)
    x = numerator / (denominator + 1e-6) - 1.0
    return 1 / (1 + np.exp(-x / slope_factor))


def _threshknot(bpp, theta=0.3, max_iter=1, allowed_buldge_len=0, min_len_helix=2):
    iteration = 0
    length = bpp.shape[0]
    bp_list = []
    new_bp = 1
    while new_bp != 0 and iteration <= max_iter:
        current_bp_list = []
        bp_list_flat = np.array(bp_list).flatten()
        if np.any(bp_list_flat):
            bpp_update = np.delete(bpp, bp_list_flat, axis=1)
            if np.any(bpp_update):
                Pmax = np.amax(bpp_update, axis=1)
        else:
            Pmax = np.amax(bpp, axis=1)
        for i in range(length):
            for j in range(i + 1, length):
                if i not in bp_list_flat and j not in bp_list_flat:
                    prob = bpp[i, j]
                    if prob == Pmax[i] and prob == Pmax[j] and prob > theta:
                        current_bp_list.append([i, j])
        new_bp = len(current_bp_list)
        iteration += 1
        if new_bp != 0 and iteration > max_iter:
            print("Reached max iteration, stopping before converged.")
        else:
            bp_list.extend(current_bp_list)

    bp_list = _check_bp_list(bp_list)
    structure = convert_bp_list_to_dotbracket(bp_list, length)
    structure = post_process_struct(structure, allowed_buldge_len, min_len_helix)
    bp_list = convert_dotbracket_to_bp_list(structure, allow_pseudoknots=True)
    return structure, bp_list


def _check_bp_list(bp_list):
    for bp in bp_list:
        bp.sort()
    bp_list.sort(key=lambda x: x[0])
    nts = [nt for bp in bp_list for nt in bp]
    if len(nts) > len(set(nts)):
        print("WARNING some nucletotides found in more than 1 bp")
        for i, bpA in enumerate(bp_list):
            for bpB in bp_list[i + 1:]:
                if bpA[0] == bpB[0] and bpA[1] == bpB[1]:
                    print("removing repeat bp", bpA)
                    bp_list = bp_list[:i] + bp_list[i + 1:]
                elif bpA[0] in bpB:
                    if abs(bpA[0] - bpA[1]) <= abs(bpB[0] - bpB[1]):
                        to_remove = bpB
                    else:
                        to_remove = bpA
                    print("WARNING base", bpA[0], "is in 2 basepairs", bpA, bpB, "THIS SHOULD BE FIXED. Removing", to_remove)
                    bp_list.remove(to_remove)
                elif bpA[1] in bpB:
                    if abs(bpA[0] - bpA[1]) <= abs(bpB[0] - bpB[1]):
                        to_remove = bpB
                    else:
                        to_remove = bpA
                    print("WARNING base", bpA[1], "is in 2 basepairs", bpA, bpB, "THIS SHOULD BE FIXED. Removing", to_remove)
                    bp_list.remove(to_remove)
    return bp_list


def _run_hotknots(seq, model="DP", param="parameters_DP03.txt"):
    hotknot_location = package_locs["hotknots"]
    cur_dir = getcwd()
    chdir(hotknot_location)
    command = [f"{hotknot_location}/HotKnots", "-noPS", "-s", seq, "-m", model, "-p", f"{hotknot_location}/params/{param}"]
    p = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = p.communicate()
    if p.returncode:
        print('ERROR: hotknots failed: on %s\n%s\n%s' % (seq, out.decode(), err.decode()))
        return ["x"*len(seq)]
    output = out.decode().split("\n")[2:-1]
    structs = []
    for struct in output:
        x = struct.split('\t')
        x2 = [x[0].split(" ")[-1], x[1]]
        structs.append(x2)
    chdir(cur_dir)
    return structs


def _ipknot_mfe(seq, model="LinearPartition-C", refinement=1, t1="auto", t2="auto"):
    """
    TODO
      -g, --gamma G             The weight for true base-pairs equivalent to 
                                '-t 1/(gamma+1)'
      -i, --allow-isolated      Allow isolated base-pairs
      -P, --param FILE          Read the energy parameter file for Vienna RNA 
                                package
      -x, --aux                 Import an auxiliary file for base-pairing 
                                probabilities
      -u, --no-levelwise        Do not perform the levelwise prediction
      -E, --energy              Output with the free energy
    """
    ipknot_location = package_locs["ipknot"]
    out_folder = get_random_folder()
    mkdir(out_folder)
    fasta_file = f"{out_folder}/temp.fasta"
    f = open(fasta_file, "w")
    f.write(">seq \n")
    f.write(seq)
    f.close()
    command = [f"{ipknot_location}/ipknot", fasta_file, "--model", model, "-r", str(refinement), "-t", str(t1), "-t", str(t2)]
    p = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = p.communicate()
    if p.returncode:
        print('ERROR: ipknot failed: on %s\n%s\n%s' % (seq, out.decode(), err.decode()))
        remove(fasta_file)
        rmdir(out_folder)
        return "x"*len(seq)
    output = out.decode().split("\n")
    remove(fasta_file)
    rmdir(out_folder)
    return output[2]


def _knotty_mfe(seq):
    knotty_location = package_locs["knotty"]
    command = [f"{knotty_location}/knotty", seq]
    p = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE,universal_newlines=True)
    try:
        out, err = p.communicate()
    except:
        print("ERROR knotty, could not communicate")
        return "x"*len(seq)
    if p.returncode:
        print('ERROR: knotty failed: on %s\n%s\n%s' % (seq, out, err))
        return "x"*len(seq)
    output = out.split("\n")
    struct = output[1].split(" ")[1]
    bp_list = convert_dotbracket_to_bp_list(struct, allow_pseudoknots=True)
    struct = convert_bp_list_to_dotbracket(bp_list, seq_len=len(struct))
    return struct


def _run_spotrna(seq, cpu=32):
    '''
    SPOT-RNA
    '''
    spotrna_location = package_locs["spotrna"]
    spotrna_conda_env = package_locs["spotrna_conda_env"]
    out_folder = get_random_folder()
    mkdir(out_folder)
    fasta_file = f"{out_folder}/temp.fasta"
    f = open(fasta_file, "w")
    f.write(">seq\n")
    f.write(seq)
    f.close()
    command = [f"{spotrna_conda_env}/python3", f"{spotrna_location}/SPOT-RNA.py", "--inputs", fasta_file, "--outputs", out_folder, "--cpu", str(cpu)]
    # keep running until output file exists
    while not path.exists(out_folder + "/seq.bpseq"):
        p = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE)
        out, err = p.communicate()
        # print(seq, out.decode(),err.decode())
        if p.returncode:
            print('ERROR: spotrna failed: on %s\n%s\n%s' % (seq, out.decode(), err.decode()))
            return "x"*len(seq)
    bp_list = bpseq_to_bp_list(out_folder + "/seq.bpseq")
    struct = convert_bp_list_to_dotbracket(bp_list, len(seq))
    bpp = prob_to_bpp(out_folder + "/seq.prob")
    remove(out_folder + "/seq.bpseq")
    remove(out_folder + "/seq.prob")
    remove(out_folder + "/seq.ct")
    remove(fasta_file)
    rmdir(out_folder)
    return struct, bpp

def _run_spotrna2(seq):
    # TODO 
    spotrna2_location = package_locs["spotrna2"]
    out_folder = get_random_folder()
    mkdir(out_folder)
    fasta_file = f"{out_folder}/temp.fasta"
    f = open(fasta_file, "w")
    f.write(">seq\n")
    f.write(seq)
    f.close()
    command = [f"{spotrna2_location}/run_spotrna2.sh", fasta_file]
    p = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = p.communicate()
    if p.returncode:
        print('ERROR: spotrna2 failed: on %s\n%s\n%s' % (seq, out.decode(), err.decode()))
        return "x"*len(seq)
    bp_list = bpseq_to_bp_list(f"{out_folder}/temp_outputs/temp.bpseq")
    struct = convert_bp_list_to_dotbracket(bp_list, len(seq))
    bpp = prob_to_bpp(f"{out_folder}/temp_outputs/temp.prob")
    for f in os.listdir(f"{out_folder}/temp_outputs"):
        remove(f)
    rmdir(f"{out_folder}/temp_outputs")
    for f in os.listdir(f"{out_folder}/temp_features"):
        remove(f)
    rmdir(f"{out_folder}/temp_features")
    remove(fasta_file)
    rmdir(out_folder)
    return struct, bpp

def _e2efold(seq):
    # only if <600
    # TODO probably plenty of options
    e2efold_location = package_locs["e2efold"]
    e2efold_conda_env = package_locs["e2efold_conda_env"]
    out_folder = get_random_folder()
    mkdir(out_folder)
    with open(f'{out_folder}/config.json', 'w') as f:
        f.write('\n'.join(['{',
                           '        "exp_name": "performance on short sequences (50-600)",',
                           f'        "test_folder": "{out_folder}/short_seqs",',
                           f'        "save_folder": "{out_folder}/short_cts",',
                           '        "gpu": "0",',
                           '        "u_net_d": 10,',
                           '        "BATCH_SIZE": 8,',
                           '        "batch_size_stage_1": 20,',
                           '        "batch_size_stage_2": 16,',
                           '        "OUT_STEP": 100,',
                           '        "LOAD_MODEL": true,',
                           '        "pp_steps": 20,',
                           '        "pp_loss": "f1",',
                           '        "pp_model": "mixed",',
                           '        "rho_per_position": "matrix",',
                           '        "data_type": "rnastralign_all_600",',
                           '        "model_type": "att_simple_fix",',
                           '        "epoches_first": 50,',
                           '        "epoches_second": 10,',
                           '        "epoches_third": 10,',
                           '        "evaluate_epi": 1,',
                           '        "evaluate_epi_stage_1": 5,',
                           '        "step_gamma": 1,',
                           '        "k": 1,',
                           '        "test": {',
                           '                "f1": true,',
                           '                "accuracy": false,',
                           '                "energy": false',
                           '        }',
                           '}']))
    mkdir(f'{out_folder}/short_seqs')
    mkdir(f'{out_folder}/short_cts')
    command = [f"{e2efold_conda_env}/python", f"{e2efold_location}/e2efold_productive_short.py", "-c", f"{out_folder}/config.json"]
    fasta_file = f"{out_folder}/short_seqs/temp.seq"
    f = open(fasta_file, "w")
    f.write(seq)
    f.close()
    # keep running until output file exists
    while not path.exists(f"{out_folder}/short_cts/temp.seq.ct"):
        out, err = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE).communicate()
    bp_list = ct_to_bp_list(f"{out_folder}/short_cts/temp.seq.ct", 1)
    struct = convert_bp_list_to_dotbracket(bp_list, len(seq))
    remove(fasta_file)
    remove(f"{out_folder}/short_cts/temp.seq.ct")
    remove(f"{out_folder}/config.json")
    rmdir(f'{out_folder}/short_seqs')
    rmdir(f'{out_folder}/short_cts')
    rmdir(out_folder)
    return struct


def _pknots(seq):
    ''' TODO
      -a          : pseudoknot approx, exclude V7-V10 and WB9-WB1
      -c          : add L^5 coaxials (V6)
      -s          : shuffle sequences
    '''
    pknots_location = package_locs["pknots"]
    out_folder = get_random_folder()
    mkdir(out_folder)
    fasta_file = f"{out_folder}/temp.fasta"
    f = open(fasta_file, "w")
    f.write(">seq \n")
    f.write(seq)
    f.close()
    outfile = f"{out_folder}/out.out"
    command = [pknots_location + "/pknots", "-k", "-g", fasta_file, outfile]
    p = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = p.communicate()
    remove(fasta_file)
    if p.returncode:
        print('ERROR: PKNOTS failed: on %s\n%s\n%s' % (seq, out.decode(), err.decode()))
        return "x"*len(seq)
    bp_list = ct_to_bp_list(outfile, 4)
    remove(outfile)
    rmdir(out_folder)
    struct = convert_bp_list_to_dotbracket(bp_list, len(seq))
    return struct


def _nupack_mfe_pk(seq):
    # TODO many nupack options... also why is this not implemented in mfe?
    nupack_location = package_locs['nupack']
    out_folder = get_random_folder()
    mkdir(out_folder)
    fasta_file = f"{out_folder}/temp"
    f = open(f'{fasta_file}.in','w')
    f.write(seq)
    f.close()
    struct = None
    command = [nupack_location+'/mfe', "-pseudo", fasta_file]
    p = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE)
    out,err = p.communicate()
    if p.returncode:
        print(f'ERROR: nupack mfe pk failed on {seq} {fasta_file} {out.decode} {err.decode}')
        return 'x'*len(seq)
    f = open(f'{fasta_file}.mfe')
    struct = f.readlines()[16][:-1]
    f.close()
    remove(f'{fasta_file}.in')
    remove(f'{fasta_file}.mfe')
    rmdir(out_folder)
    return struct

