from arnie.mfe_bootstrap import mfe_bootstrap
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys,ast

def evaluate_shapeknots_v_fold(seq,reactivity,chemical_probe="shape",num_bootstrap=100,Pcutoff=0.3):

    if chemical_probe=="shape":
        _, shapeknots_bpp1 = mfe_bootstrap(seq, num_bootstrap, package='rnastructure', T=37,
            shape_signal=reactivity, pk=True)
        _, shapeknots_bpp2 = mfe_bootstrap(seq, num_bootstrap, package='rnastructure', T=37,
            shape_signal=reactivity, pk=True)
        _, shapeknots_bpp3 = mfe_bootstrap(seq, num_bootstrap, package='rnastructure', T=37,
            shape_signal=reactivity, pk=True)
        _, shapeknots_bpp4 = mfe_bootstrap(seq, num_bootstrap, package='rnastructure', T=37,
            shape_signal=reactivity, pk=True)
        _, shapeknots_bpp5 = mfe_bootstrap(seq, num_bootstrap, package='rnastructure', T=37,
            shape_signal=reactivity, pk=True)
        _, shapeknots_bpp6 = mfe_bootstrap(seq, num_bootstrap, package='rnastructure', T=37,
            shape_signal=reactivity, pk=True)
        _, fold_bpp1 = mfe_bootstrap(seq, num_bootstrap, package='rnastructure', T=37,
            shape_signal=reactivity, pk=False)
        _, fold_bpp2 = mfe_bootstrap(seq, num_bootstrap, package='rnastructure', T=37,
            shape_signal=reactivity, pk=False)
        _, fold_bpp3 = mfe_bootstrap(seq, num_bootstrap, package='rnastructure', T=37,
            shape_signal=reactivity, pk=False)
        np.save(seq+"_bpp_shapeknot_1.npy", shapeknots_bpp1)
        np.save(seq+"_bpp_shapeknot_2.npy", shapeknots_bpp2)
        np.save(seq+"_bpp_shapeknot_3.npy", shapeknots_bpp3)
        np.save(seq+"_bpp_shapeknot_4.npy", shapeknots_bpp4)
        np.save(seq+"_bpp_shapeknot_5.npy", shapeknots_bpp5)
        np.save(seq+"_bpp_shapeknot_6.npy", shapeknots_bpp6)
        np.save(seq+"_bpp_fold_1.npy", fold_bpp1)
        np.save(seq+"_bpp_fold_2.npy", fold_bpp2)
        np.save(seq+"_bpp_fold_3.npy", fold_bpp3)
        result = evaluate_per_bp_pairwise_mean_diff(seq,[shapeknots_bpp1,shapeknots_bpp2,shapeknots_bpp3,
            shapeknots_bpp4,shapeknots_bpp5,shapeknots_bpp6],[fold_bpp1,fold_bpp2,fold_bpp3],Pcutoff)

    elif chemical_probe=="dms":
        _, shapeknots_bpp1 = mfe_bootstrap(seq, num_bootstrap, package='rnastructure', T=37,
            dms_signal=reactivity, pk=True)
        _, shapeknots_bpp2 = mfe_bootstrap(seq, num_bootstrap, package='rnastructure', T=37,
            dms_signal=reactivity, pk=True)
        _, shapeknots_bpp3 = mfe_bootstrap(seq, num_bootstrap, package='rnastructure', T=37,
            dms_signal=reactivity, pk=True)
        _, shapeknots_bpp4 = mfe_bootstrap(seq, num_bootstrap, package='rnastructure', T=37,
            dms_signal=reactivity, pk=True)
        _, shapeknots_bpp5 = mfe_bootstrap(seq, num_bootstrap, package='rnastructure', T=37,
            dms_signal=reactivity, pk=True)
        _, shapeknots_bpp6 = mfe_bootstrap(seq, num_bootstrap, package='rnastructure', T=37,
            dms_signal=reactivity, pk=True)
        _, fold_bpp1 = mfe_bootstrap(seq, num_bootstrap, package='rnastructure', T=37,
            dms_signal=reactivity, pk=False)
        _, fold_bpp2 = mfe_bootstrap(seq, num_bootstrap, package='rnastructure', T=37,
            dms_signal=reactivity, pk=False)
        _, fold_bpp3 = mfe_bootstrap(seq, num_bootstrap, package='rnastructure', T=37,
            dms_signal=reactivity, pk=False)
        np.save(seq+"_bpp_shapeknot_1.npy", shapeknots_bpp1)
        np.save(seq+"_bpp_shapeknot_2.npy", shapeknots_bpp2)
        np.save(seq+"_bpp_shapeknot_3.npy", shapeknots_bpp3)
        np.save(seq+"_bpp_shapeknot_4.npy", shapeknots_bpp4)
        np.save(seq+"_bpp_shapeknot_5.npy", shapeknots_bpp5)
        np.save(seq+"_bpp_shapeknot_6.npy", shapeknots_bpp6)
        np.save(seq+"_bpp_fold_1.npy", fold_bpp1)
        np.save(seq+"_bpp_fold_2.npy", fold_bpp2)
        np.save(seq+"_bpp_fold_3.npy", fold_bpp3)
        result = evaluate_per_bp_pairwise_mean_diff(seq,[shapeknots_bpp1,shapeknots_bpp2,shapeknots_bpp3,
            shapeknots_bpp4,shapeknots_bpp5,shapeknots_bpp6],[fold_bpp1,fold_bpp2,fold_bpp3],Pcutoff)

    else:
        print("ERROR: chemical probe not recognize, must be shape or dms")
    
    return pd.DataFrame(result,columns=["seq","bp","inter/intra","diff"])

    
def evaluate_per_bp_pairwise_mean_diff(seq,bppAs_,bppBs_,Pcutoff):
    bppAs = bppAs_.copy()
    bppBs = bppBs_.copy()
    compare_locations = np.argwhere((bppAs[0]>Pcutoff) | (bppAs[1]>Pcutoff) | (bppAs[2]>Pcutoff) |
                                    (bppAs[3]>Pcutoff) | (bppAs[4]>Pcutoff) | (bppAs[5]>Pcutoff) |
                                    (bppBs[0]>Pcutoff) | (bppBs[1]>Pcutoff) | (bppBs[2]>Pcutoff))
    cutoff_filter = ((bppAs[0]>Pcutoff) | (bppAs[1]>Pcutoff) | (bppAs[2]>Pcutoff) |
                                    (bppAs[3]>Pcutoff) | (bppAs[4]>Pcutoff) | (bppAs[5]>Pcutoff) |
                                    (bppBs[0]>Pcutoff) | (bppBs[1]>Pcutoff) | (bppBs[2]>Pcutoff))
    for i in range(3):
        bppBs[i][~cutoff_filter] = 0
    for i in range(6):
        bppAs[i][~cutoff_filter] = 0
    
    diff_intra = []
    diff_inter = []

    for i in range(5):
        for j in range(i+1,6):
            diff_intra.append(abs(bppAs[i]-bppAs[j]))
    for i in range(6):
        for j in range(3):
            diff_inter.append(abs(bppAs[i]-bppBs[j]))
    
    results = []
    
    for loc in compare_locations:
        for diff in diff_intra:
            results.append([seq,str(loc),"intra",diff[loc[0],loc[1]]])
        for diff in diff_inter:
            results.append([seq,str(loc),"inter",diff[loc[0],loc[1]]])
    return results

if __name__=='__main__':
    seq = sys.argv[1]
    react_file = sys.argv[2]
    with open(react_file) as f:
        reactivity = [float(line) for line in f]
    data = evaluate_shapeknots_v_fold(seq,reactivity,chemical_probe="shape",num_bootstrap=10,Pcutoff=0.3)
    sns.violinplot(x="bp", y="diff",
                hue="inter/intra", data=data[data.seq==seq])
    plt.savefig(seq+".png")
