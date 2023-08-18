from arnie.pk_predictors import pk_predict, pk_predict_from_bpp
from arnie.utils import prob_to_bpp, load_package_locations
import numpy as np

# TODO e2efold is stochastic?
# TODO spotrna2 add in?

samiv_seq = "GGUCAUGAGUGCCAGCGUCAAGCCCCGGCUUGCUGGCCGGCAACCCUCCAACCGCGGUGGGGUGCCCCGGGUGAUGACCAGGUUGAGUAGCCGUGACGGCUACGCGGCAAGCGCGGGUC"
samiv_struct = "((((....(.((((((....((.[[[[[)).)))))))(((..((((((..{{{{)).)))).)))]]]]]....))))..((((.(((((.......))))).))))....}}}}..."

pk_res = {"hotknots": "((((((....(((((((......[[[[[..))))))).....(((((....{{{{{..)))))...]]]]]..))))))..((((.(((((((...))))))).))))...}}}}}...",
          "ipknot": "[[[[[[....[[[..........(((((((((((]]].((((.(((.((......)).))).)))).......]]]]]].......(((((((...)))))))..))))))).))))..",
          "knotty": "(((..[[[[[))).]].]]]...((((((((((((...((((.(((.((......)).))).))))((.(((....))).))....(((((((...))))))).)))))))).))))..",
          "spotrna": "(((((((....(((((...((.....[[.)))))))]](((.((((((((....)))).)))))))......))))))).......((((((((.))))))))(.........).....",
          "e2efold": ".....(............((...(.......(.............)..........)).)........(.(......(.(.(.([...).]).)...)....).).....)........",
          "pknots": "(((.......))).((.....))(((((((((((((((((((.((((((......)).))))))))...(((....))).)))...((((((.....)))))).)))))))).)))).."}

# threshknot_theta_maxIter_buldge_helix
# hungarian_theta_buldge_helix_exp_sig_0p_1p_ln_unpaired
bpp_heuristics = {"threshknot_0.1_1_0_1": "...((......(...........[[[[[[[[[[[[)..((((.(((.((....())).))).))))[[.[[[..))]]].]]....(((((((...))))))).]]]]]]]].]]]]..",
                  "threshknot_0.4_1_0_1": "...........................(((((((....(((......(........)......)))....................(((((((...)))))))..))))))).......",
                  "threshknot_0.9_1_0_1": ".......................................................................................................................",
                  "threshknot_0.1_1_0_3": ".......................((((((((((((...((((.(((............))).))))...(((....))).......(((((((...))))))).)))))))).))))..",
                  "threshknot_0.1_5_0_1": "(..[[.....[[)........(.((((((((((((]].((((.(((.((....())).))).))))((.(((..]]))).))....(((((((...))))))).)))))))).)))).)",
                  "hungarian_0.3_0_1_1_None_0.1_0.9_False_True": ".......................(((((((((((....(((..(((.((......)).)))..)))....................(((((((...)))))))..))))))).))))..",
                  "hungarian_0.3_0_2_4_None_0.1_0.9_False_True": ".......................(((((((((((....((((.(((.((......)).))).))))....................(((((((...)))))))..))))))).))))..",
                  "hungarian_0.3_0_1_1_None_0_1_False_False": ".......................(((((((((((....((((.(((.((......)).))).))))....................(((((((...)))))))..))))))).))))..",
                  "hungarian_0.3_0_1_1_3_0.1_0.9_False_True": ".......................(((((((((((....(((..(((.((......)).)))..)))....................(((((((...)))))))..))))))).))))..",
                  "hungarian_0.8_0_1_1_None_0.1_0.9_False_True": "......................................................................................((((((.....))))))................"}


def test_pk(pkg):
    print("Testing", pkg)
    pred = pk_predict(samiv_seq, pkg)

    assert(pred == pk_res[pkg])


# def bpps and output expected
bpp_file = "test_files/samiv_eternafold.prob"
bpp = prob_to_bpp(bpp_file)


def test_pk_from_bpp():
    print("Testing threshknot")
    assert(bpp_heuristics["threshknot_0.1_1_0_1"] == pk_predict_from_bpp(bpp, heuristic="threshknot", theta=0.1, max_iter=1, allowed_buldge_len=0, min_len_helix=1))
    assert(bpp_heuristics["threshknot_0.4_1_0_1"] == pk_predict_from_bpp(bpp, heuristic="threshknot", theta=0.4, max_iter=1, allowed_buldge_len=0, min_len_helix=1))
    assert(bpp_heuristics["threshknot_0.9_1_0_1"] == pk_predict_from_bpp(bpp, heuristic="threshknot", theta=0.9, max_iter=1, allowed_buldge_len=0, min_len_helix=1))
    assert(bpp_heuristics["threshknot_0.1_1_0_3"] == pk_predict_from_bpp(bpp, heuristic="threshknot", theta=0.1, max_iter=1, allowed_buldge_len=0, min_len_helix=3))
    assert(bpp_heuristics["threshknot_0.1_5_0_1"] == pk_predict_from_bpp(bpp, heuristic="threshknot", theta=0.1, max_iter=5, allowed_buldge_len=0, min_len_helix=1))
    print("Testing hungarian")
    assert(bpp_heuristics["hungarian_0.3_0_1_1_None_0.1_0.9_False_True"] == pk_predict_from_bpp(bpp, heuristic="hungarian", theta=0.3, allowed_buldge_len=0, min_len_helix=1,
                                                                                                exp=1, sigmoid_slope_factor=None, prob_to_0_threshold_prior=0.1, prob_to_1_threshold_prior=0.9, ln=False, add_p_unpaired=True))
    assert(bpp_heuristics["hungarian_0.3_0_2_4_None_0.1_0.9_False_True"] == pk_predict_from_bpp(bpp, heuristic="hungarian", theta=0.3, allowed_buldge_len=2, min_len_helix=4,
                                                                                                exp=1, sigmoid_slope_factor=None, prob_to_0_threshold_prior=0.1, prob_to_1_threshold_prior=0.9, ln=False, add_p_unpaired=True))
    assert(bpp_heuristics["hungarian_0.3_0_1_1_None_0_1_False_False"] == pk_predict_from_bpp(bpp, heuristic="hungarian", theta=0.3, allowed_buldge_len=0, min_len_helix=1,
                                                                                             exp=1, sigmoid_slope_factor=None, prob_to_0_threshold_prior=0, prob_to_1_threshold_prior=1, ln=False, add_p_unpaired=False))
    assert(bpp_heuristics["hungarian_0.3_0_1_1_3_0.1_0.9_False_True"] == pk_predict_from_bpp(bpp, heuristic="hungarian", theta=0.3, allowed_buldge_len=0, min_len_helix=1,
                                                                                             exp=1, sigmoid_slope_factor=3, prob_to_0_threshold_prior=0.1, prob_to_1_threshold_prior=0.9, ln=False, add_p_unpaired=True))
    assert(bpp_heuristics["hungarian_0.8_0_1_1_None_0.1_0.9_False_True"] == pk_predict_from_bpp(bpp, heuristic="hungarian", theta=0.8, allowed_buldge_len=0, min_len_helix=1,
                                                                                                exp=1, sigmoid_slope_factor=None, prob_to_0_threshold_prior=0.1, prob_to_1_threshold_prior=0.9, ln=False, add_p_unpaired=True))


if __name__ == '__main__':
    package_locs = load_package_locations()
    pk_predictors = ["spotrna", "e2efold", "hotknots", "ipknot", "knotty", "pknots"]
    for pkg in pk_predictors:
        if pkg not in package_locs:
            print("Warning:", pkg, "is not found in the ARNIEFILE, not testing.")
        else:
            test_pk(pkg)
    test_pk_from_bpp()
