import pandas as pd
from arnie.utils import is_PK, compare_structure_to_native, compare_structures_to_natives, post_process_struct
from arnie.pk_predictors import pk_predict, pk_predict_from_bpp
from arnie.bpps import bpps

# example
# evaluate_pk_predictors('test.csv',pk_predictors=["hotknots","ipknot"],helix_filters=[[4,1]],bpp_packages=["contrafold"],
#		from_bpp_heuristic=["threshknot","hungarian"],from_bpp_heuristic_options={"threshknot":[{'theta':0.2},{'theta':0.3},{'theta':0.9}]})


def evaluate_pk_predictors(csv_file, pk_predictors=[], pk_predictors_options={}, helix_filters=[], bpp_packages=[], bpp_packages_options={},
                           bpp_linear_packages=[], bpp_linear_packages_options=[], from_bpp_heuristic=[], from_bpp_heuristic_options={}):
    # structure, sequence columns
    # pk_predictors: list, allowed: ["hotknots", "ipknot", "knotty", "spotrna", "e2efold", "pknots"]
    # helix_filters: [[h1,b1],[h2,b2]] where all helices under h in length are delted allowing for b bugles (eg if b=2 allows 0-1,1-1,2-1,and 2-2 buldges)
    # bpp_packages: list, allowed: ['nupack','rnastructure','rnasoft','eternafold','vienna','contrafold']
    # bpp_linear_packages: list allowed: ['eternafold','vienna','contrafold']
    # from_bpp_heuristic: list allowed, ['threshknot','hungarian']
    # the options arguments are dictionary for each respective predictor that has a list of dictionaries for the various option sets
    # eg if want to run threshknot with different threshold: from_bpp_heuristic_options={"threshknot":[{'theta':0.2},{'theta':0.3},{'theta':0.9}]
    data = pd.read_csv(csv_file)

    # run predictors
    for predictor in pk_predictors:
        if predictor in pk_predictors_options:
            for option_set in pk_predictors_options[predictor]:
                data[predictor + "_" + "_".join([str(x) + str(y) for x, y in option_set.items()]) + "_structure"] = data.apply(lambda x: pk_predict(x["sequence"], predictor, **option_set), axis=1)
        else:
            data[predictor + "_structure"] = data.apply(lambda x: pk_predict(x["sequence"], predictor), axis=1)

    # run bpp predictors
    all_bpps = {}
    for pkg in bpp_packages:
        if pkg in bpp_packages_options:
            for option_set in bpp_packages_options[pkg]:
                data[pkg + "_" + "_".join([str(x) + str(y) for x, y in option_set.items()]) + "_bpp"] = data.apply(lambda x: bpps(x["sequence"], package=pkg, **option_set), axis=1)
        else:
            data[pkg + "_bpp"] = data.apply(lambda x: bpps(x["sequence"], package=pkg), axis=1)
    for pkg in bpp_linear_packages:
        if pkg in bpp_linear_packages_options:
            for option_set in bpp_linear_packages_options[pkg]:
                data[pkg + "_" + "_".join([str(x) + str(y) for x, y in option_set.items()]) + "_linear_bpp"] = data.apply(lambda x: bpps(x["sequence"], linear=True, package=pkg, **option_set), axis=1)
        else:
            data[pkg + "_linear_bpp"] = data.apply(lambda x: bpps(x["sequence"], package=pkg, linear=True), axis=1)
    bpp_cols = [x for x in list(data.columns) if "_bpp" in x]
    for heuristic in from_bpp_heuristic:
        if heuristic in from_bpp_heuristic_options:
            for option_set in from_bpp_heuristic_options[heuristic]:
                for bpp_col in bpp_cols:
                    data[bpp_col + "_" + heuristic + "_" + "_".join([str(x) + str(y) for x, y in option_set.items()]) + "_structure"] = data.apply(lambda x: pk_predict_from_bpp(x[bpp_col], heuristic=heuristic, **option_set), axis=1)
        else:
            for bpp_col in bpp_cols:
                data[bpp_col + "_" + heuristic + "_structure"] = data.apply(lambda x: pk_predict_from_bpp(x[bpp_col], heuristic=heuristic), axis=1)
    # save bpp in readable format once csv written
    # to read use: np.array(ast.literal_eval(row.bpp))
    for bpp_col in bpp_cols:
        data[bpp_col] = data.apply(lambda x: [bpp.tolist() for bpp in x[bpp_col]], axis=1)

    # post process filters
    structure_cols = [x for x in list(data.columns) if "structure" in x]
    for structure_col in structure_cols:
        for h, b in helix_filters:
            data[structure_col + "_filter_<" + str(h) + "helix_allow" + str(b) + "buldge"] = data.apply(lambda x: post_process_struct(x[structure_col], allowed_buldge_len=b, min_len_helix=h), axis=1)

    # for every structure add column for whether it is a PK
    structure_cols = [x for x in list(data.columns) if "structure" in x]

    for structure_col in structure_cols:
        data[structure_col + "_is_PK"] = data.apply(lambda x: is_PK(x[structure_col]), axis=1)

    # for every structure add columns for sensitivity, ppv, and f1 score (pk_only, no_pk, all baspairs)
    for structure_col in structure_cols:
        data[structure_col + "_sensitivity"] = data.apply(lambda x: compare_structure_to_native(x[structure_col], x["structure"], metric="sensitivity"), axis=1)
        data[structure_col + "_PPV"] = data.apply(lambda x: compare_structure_to_native(x[structure_col], x["structure"], metric="PPV"), axis=1)
        data[structure_col + "_F1_score"] = data.apply(lambda x: compare_structure_to_native(x[structure_col], x["structure"], metric="F1_score"), axis=1)
        data[structure_col + "_PKonly_sensitivity"] = data.apply(lambda x: compare_structure_to_native(x[structure_col], x["structure"], metric="sensitivity", PK_involved=True), axis=1)
        data[structure_col + "_PKonly_PPV"] = data.apply(lambda x: compare_structure_to_native(x[structure_col], x["structure"], metric="PPV", PK_involved=True), axis=1)
        data[structure_col + "_PKonly_F1_score"] = data.apply(lambda x: compare_structure_to_native(x[structure_col], x["structure"], metric="F1_score", PK_involved=True), axis=1)
        data[structure_col + "_noPK_sensitivity"] = data.apply(lambda x: compare_structure_to_native(x[structure_col], x["structure"], metric="sensitivity", PK_involved=False), axis=1)
        data[structure_col + "_noPK_PPV"] = data.apply(lambda x: compare_structure_to_native(x[structure_col], x["structure"], metric="PPV", PK_involved=False), axis=1)
        data[structure_col + "_noPK_F1_score"] = data.apply(lambda x: compare_structure_to_native(x[structure_col], x["structure"], metric="F1_score", PK_involved=False), axis=1)

    # output full sensitivty, ppv, f1 for al, binary
    natives = data["structure"].tolist()
    overview_results = {}
    for structure_col in structure_cols:
        overview_results[structure_col] = {}
        structs = data[structure_col].tolist()
        overview_results[structure_col]["all_bps"] = compare_structures_to_natives(structs, natives, comparison="basepairs", metric="all")
        overview_results[structure_col]["PKonly_bps"] = compare_structures_to_natives(structs, natives, comparison="PK_basepairs", metric="all")
        overview_results[structure_col]["noPK_bps"] = compare_structures_to_natives(structs, natives, comparison="non_PK_basepairs", metric="all")
        overview_results[structure_col]["is_PK"] = compare_structures_to_natives(structs, natives, comparison="is_PK", metric="all")

    # save csv
    data.to_csv(".".join(csv_file.split(".")[:-1]) + "_output.csv", index=False)
    overview_results = pd.DataFrame.from_dict({(i, j): overview_results[i][j] for i in overview_results.keys() for j in overview_results[i].keys()}, orient='index')
    overview_results.to_csv(".".join(csv_file.split(".")[:-1]) + "_summary.csv", index_label=['predictor', 'evaluation'])
