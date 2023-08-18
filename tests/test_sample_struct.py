from arnie.sample_structures import sample_structures

sample_seq = 'GGGGAAAACCCC'


def test_sample_seq():

    struct_list = sample_structures(
        sample_seq, n_samples=10, package='vienna_2')
    # sample structures no longer returns energy or prob?
    # print(ener_list) # , ener_list, prob_list
    # print(prob_list)
    return


if __name__ == '__main__':
    test_sample_seq()
    # test_pkg_w_bpps(pkg.lower())
