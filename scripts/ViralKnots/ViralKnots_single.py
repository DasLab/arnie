from ViralKnots_utils import *

if __name__=='__main__':

    parser=argparse.ArgumentParser()

    parser.add_argument('--seqs', nargs='+', required=True)
    parser.add_argument('--coords', nargs='+', required=True)
    parser.add_argument('--window', '-w', type=int, required=True)
    parser.add_argument('--pk_predictors', nargs='+', required=True)
    parser.add_argument('--bpp_package', default='contrafold', type=str)
    parser.add_argument('--temp_folder', default=None, type=str)
    parser.add_argument('--linear_partition', default=None, type=bool)
    args=parser.parse_args()

    struct_list = []
    for pk_predictor in args.pk_predictors:
        for seq,coord in zip(args.seqs, args.coords):
            struct_list.append(get_structure(seq, coord, pk_predictor, args.window, args.bpp_package, args.linear_partition))

    df = pd.DataFrame(struct_list,columns=["predictor","start","end","sequence", "struct", "pseudoknot"])

    csv_name = ''
    for i in range(len(args.coords)):
        csv_name += args.coords[i] + '_'

    df.to_csv(f"{args.temp_folder}/{csv_name}.csv", index=False)

