from arnie.utils import *


s = "(((....)).)...(..)....(((..(((....))))))"
s_0_2 = ".((....)).............(((..(((....))))))"
s_0_3 = "......................(((..(((....))))))"
s_1_3 = "(((....)).)...........(((..(((....))))))"
s_2_3 = "(((....)).)...........(((..(((....))))))"
s_0_4 = "........................................"
s_1_4 = "........................................"
s_2_4 = "......................(((..(((....))))))"
s_1_2 = "(((....)).)...........(((..(((....))))))"

s_all_helices = [[[0, 10]],
                 [[1, 8], [2, 7]],
                 [[14, 17]],
                 [[22, 39], [23, 38], [24, 37]],
                 [[27, 36], [28, 35], [29, 34]]]
s_1_helices = [[[0, 10], [1, 8], [2, 7]],
               [[14, 17]],
               [[22, 39], [23, 38], [24, 37]],
               [[27, 36], [28, 35], [29, 34]]]
s_2_helices = [[[0, 10], [1, 8], [2, 7]],
               [[14, 17]],
               [[22, 39], [23, 38], [24, 37], [27, 36], [28, 35], [29, 34]]]


pk = "(((.((([..[[..))))((...)){...]]]...)})"
pk_0_2 = "....(((...[[..))).((...))....]]......."
pk_0_3 = "....(((.......)))....................."
pk_1_3 = "..(.(((.......))))...................."
pk_2_3 = "..(.((([..[[..))))...........]]]......"
pk_0_4 = "......................................"
pk_1_4 = "..(.(((.......))))...................."
pk_2_4 = "..(.(((.......))))...................."
pk_1_2 = "(((.(((...[[..))))((...))....]]....).)"

pk_all_helices = [[[0, 37]],
                  [[1, 35]],
                  [[2, 17]],
                  [[4, 16], [5, 15], [6, 14]],
                  [[7, 31]],
                  [[10, 30], [11, 29]],
                  [[18, 24], [19, 23]],
                  [[25, 36]]]
pk_1_helices = [[[0, 37], [1, 35]],
                [[2, 17], [4, 16], [5, 15], [6, 14]],
                [[7, 31]],
                [[10, 30], [11, 29]],
                [[18, 24], [19, 23]],
                [[25, 36]]]
pk_2_helices = [[[0, 37], [1, 35]],
                [[2, 17], [4, 16], [5, 15], [6, 14]],
                [[7, 31], [10, 30], [11, 29]],
                [[18, 24], [19, 23]],
                [[25, 36]]]


def test_getting_helix():
    assert(get_helices(s, allowed_buldge_len=0) == s_all_helices)
    assert(get_helices(pk, allowed_buldge_len=0) == pk_all_helices)
    assert(get_helices(s, allowed_buldge_len=1) == s_1_helices)
    assert(get_helices(pk, allowed_buldge_len=1) == pk_1_helices)
    assert(get_helices(s, allowed_buldge_len=2) == s_2_helices)
    assert(get_helices(pk, allowed_buldge_len=2) == pk_2_helices)


def test_removing_helix():
    assert(post_process_struct(s, allowed_buldge_len=0, min_len_helix=1) == s)
    # note PKs may swap around their bracket types so fairest to compare bp_list always!
    assert(convert_dotbracket_to_bp_list(post_process_struct(pk, allowed_buldge_len=0, min_len_helix=1), len(pk)) == convert_dotbracket_to_bp_list(pk, len(pk)))
    assert(post_process_struct(s, allowed_buldge_len=0, min_len_helix=2) == s_0_2)
    assert(post_process_struct(pk, allowed_buldge_len=0, min_len_helix=2) == pk_0_2)
    assert(post_process_struct(s, allowed_buldge_len=0, min_len_helix=3) == s_0_3)
    assert(post_process_struct(pk, allowed_buldge_len=0, min_len_helix=3) == pk_0_3)
    assert(post_process_struct(s, allowed_buldge_len=1, min_len_helix=3) == s_1_3)
    assert(post_process_struct(pk, allowed_buldge_len=1, min_len_helix=3) == pk_1_3)
    assert(post_process_struct(s, allowed_buldge_len=2, min_len_helix=3) == s_2_3)
    assert(post_process_struct(pk, allowed_buldge_len=2, min_len_helix=3) == pk_2_3)
    assert(post_process_struct(s, allowed_buldge_len=0, min_len_helix=4) == s_0_4)
    assert(post_process_struct(pk, allowed_buldge_len=0, min_len_helix=4) == pk_0_4)
    assert(post_process_struct(s, allowed_buldge_len=1, min_len_helix=4) == s_1_4)
    assert(post_process_struct(pk, allowed_buldge_len=1, min_len_helix=4) == pk_1_4)
    assert(post_process_struct(s, allowed_buldge_len=2, min_len_helix=4) == s_2_4)
    assert(post_process_struct(pk, allowed_buldge_len=2, min_len_helix=4) == pk_2_4)
    assert(post_process_struct(pk, allowed_buldge_len=1, min_len_helix=2) == pk_1_2)
    assert(post_process_struct(s, allowed_buldge_len=1, min_len_helix=2) == s_1_2)

if __name__ == '__main__':
    test_getting_helix()
    test_removing_helix()
