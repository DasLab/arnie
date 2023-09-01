from arnie.utils import *

samiv_struct = "((((....(.((((((....((.[[[[[)).)))))))(((..((((((..{{{{)).)))).)))]]]]]....))))..((((.(((((.......))))).))))....}}}}..."

hotknots = "((((((....(((((((......[[[[[..))))))).....(((((....{{{{{..)))))...]]]]]..))))))..((((.(((((((...))))))).))))...}}}}}..."
ipknot = "[[[[[[....[[[..........(((((((((((]]].((((.(((.((......)).))).)))).......]]]]]].......(((((((...)))))))..))))))).)))).."
knotty = "(((..[[[[[))).]].]]]...[[[[[[[[[[[[...[[[[.[[[.[[......]].]]].]]]]((.(((....))).))....(((((((...))))))).]]]]]]]].]]]].."
spotrna = ".(((((((....(((((...((.....[[.)))))))]](((.((((((((....)))).)))))))......))))))).......((((((((.))))))))(.........)...."
e2efold = "......(............((...(.......(.............)..........)).)........(.(......(.(.(.([...).]).)...)....).).....)......."
pknots = ".(((.......))).((.....))(((((((((((((((((((.((((((......)).))))))))...(((....))).)))...((((((.....)))))).)))))))).))))."
empty = "." * len(samiv_struct)


def test_is_pk():
    assert(is_PK(samiv_struct))
    assert(not is_PK(pknots))


def test_compare_struct():
    assert(compare_structure_to_native(hotknots, samiv_struct, metric="PPV") == 0.8205128205128205)
    assert(compare_structure_to_native(hotknots, samiv_struct, metric="sensitivity") == 0.8)
    assert(compare_structure_to_native(hotknots, samiv_struct, metric="F1_score") == 0.810126582278481)
    assert(compare_structure_to_native(hotknots, samiv_struct, metric="all")["F1_score"] == 0.810126582278481)
    assert(compare_structure_to_native(empty, samiv_struct, metric="all")["F1_score"] == 0)
    assert(compare_structure_to_native(hotknots, samiv_struct, metric="all", PK_involved=True)["F1_score"] == 0.7796610169491526)
    assert(compare_structure_to_native(hotknots, samiv_struct, metric="all", PK_involved=False)["F1_score"] == 0.9)


def test_compare_structs():
    assert(0.4266666666666667 == compare_structures_to_natives([hotknots, spotrna], [samiv_struct, samiv_struct], comparison="basepairs")['PPV'])
    assert(1.0 == compare_structures_to_natives([hotknots, spotrna], [samiv_struct, samiv_struct], comparison="is_PK")["F1_score"])
    assert(1.0 == compare_structures_to_natives([hotknots, spotrna, pknots, empty], [samiv_struct, samiv_struct, samiv_struct, samiv_struct], comparison="is_PK", metric="PPV"))
    assert(0.5 == compare_structures_to_natives([hotknots, spotrna, pknots, empty], [samiv_struct, samiv_struct, samiv_struct, samiv_struct], comparison="is_PK", metric="sensitivity"))
    assert(0.6666666666666666 == compare_structures_to_natives([hotknots, spotrna, pknots, empty], [samiv_struct, samiv_struct, samiv_struct, samiv_struct], comparison="is_PK", metric="F1_score"))
    assert(0.25 == compare_structures_to_natives([hotknots, spotrna, pknots, empty], [samiv_struct, samiv_struct, samiv_struct, samiv_struct], comparison="non_PK_basepairs")["sensitivity"])
    assert(compare_structures_to_natives([hotknots, spotrna, pknots, empty], [samiv_struct, samiv_struct, samiv_struct, samiv_struct], comparison="PK_basepairs")["F1_score"] == 0.28571428571428575)


if __name__ == '__main__':
    test_is_pk()
    test_compare_struct()
    test_compare_structs()
