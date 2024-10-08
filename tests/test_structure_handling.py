import pytest
from arnie import utils

success_cases = [
  ["......", False],
  ["......", True],
  ["(((((......)))))", False],
  ["(((((......)))))", True],
  ["[[[[[......]]]]]", False],
  ["{{{{{......}}}}}", False],
  ["<<<<<......>>>>>", False],
  ["((((((((...........)).))))))", False],
  ["((((((((...........)).))))))", True],
  ["(((.((((((((..((((.(((((....)).)))..))))..))))((((...))))))))...)))", False],
  ["(((.((((((((..((((.(((((....)).)))..))))..))))((((...))))))))...)))", True],
  [".(((((((((.((((....))))....((((....))))..))))..(((((......)))))...(((.((.(((((((((((((...((((..((((.....))))...)))).))))))))))))))).)))..)))))..", False],
  [".(((((((((.((((....))))....((((....))))..))))..(((((......)))))...(((.((.(((((((((((((...((((..((((.....))))...)))).))))))))))))))).)))..)))))..", True],
  [".....((((((.....))))))....(((((([[[[[[[[[[[........))))))]]]]]]]]]]]........(((((((....))))))).....................", True],
  ["(((.[[[.(((...))).]]].)))", True],
  ["(((..[[[.(((...))))))]]]", True],
  ["([{<a.aaa....)]}>AAA.A", True],
  [".....[[[[[[.....]]]]]]....[[[[[[[[[[[....]]]]]]]]]]].......(((((((((((.[[[[[)))))))))))[[[[[[[[[[[[....]]]]].]]]]]]]....]]]]].[[[[[[[[[[[[[....]]]]]]]]]]]]].....................", True],
  [".....((((((.....))))))....(((((((((((....)))))))))))........((((((((((.<<<<<))))))))))(((((((((((((....))))).))))))))...>>>>>.(((((((((((((....))))))))))))).....................", True],
  ["((.(.(.((.(.....).)).)....(((((((((((....)))))))))))......)((((((((((...{{{{.)))))))))).(((((((((((....))))).))))))...)[}}}}..[[[[[[[[[[[[[)...]]]]]]]]]]]]]...].................", True],
  ["(((((((([.{)][..((.{]).)(.}).})))))))", True],
  ["-(((((..(((((....)))))((((((((.((...-----((((((..(((((((..)))))))(((((({..[[[[[[[)))))))))))..).))).))))))))))))}.]]]]]]]...", True],
  ["((((({<A[[[....))))).......}>]a]]", True],
  [".....(.((.(.....).)).)....(.((..({(<.{.a.a{........(..(((((.(.[....).)))).)...)).).}).((.(}...}...)>)).A.A...).)..............(((((((((((.(....).)))))))))))............]........", True],
  ["..((...((.........))......(((((((....))))))).(((((((....))))))).(((((((....))))))).(((((((((....)))))))))..(((((((....)))))))..(((((((....)))))))...[[[[[[[......[[[[[[[[......[[.....)).....((..........]]......]]]]]]]]......................))..............]]]]]]].....(((((((....)))))))......................", True],
  ["((((((((((((((((....)))))))((((.(((((((((((......)))))[[[[[[.)))))).))))(((((((((((((.((.(.....(......(((.(((((((((((.((((.(((((..]]]]]].))))).)))).))(((((((....)))))))(((((((((....)))))))(((((((....)))))))(((.(.((((((((((......)))))[[[[[[.))))).).)))))))))))))).)))..............)..).)).)))))))))))(((((((.(((((..]]]]]].))))).)))))))(((((((....))))))))))))))))))", True],
  ["(((....))) (((....)))", False],
  ["(((....))) (((....)))", True],
]

success_expected_output = [
  "......",
  "......",
  "(((((......)))))",
  "(((((......)))))",
  "(((((......)))))",
  "(((((......)))))",
  "(((((......)))))",
  "((((((((...........)).))))))",
  "((((((((...........)).))))))",
  "(((.((((((((..((((.(((((....)).)))..))))..))))((((...))))))))...)))",
  "(((.((((((((..((((.(((((....)).)))..))))..))))((((...))))))))...)))",
  ".(((((((((.((((....))))....((((....))))..))))..(((((......)))))...(((.((.(((((((((((((...((((..((((.....))))...)))).))))))))))))))).)))..)))))..",
  ".(((((((((.((((....))))....((((....))))..))))..(((((......)))))...(((.((.(((((((((((((...((((..((((.....))))...)))).))))))))))))))).)))..)))))..",
  ".....((((((.....))))))....(((((([[[[[[[[[[[........))))))]]]]]]]]]]]........(((((((....))))))).....................",
  "(((.(((.(((...))).))).)))",
  "(((..[[[.(((...))))))]]]",
  "([{<a.aaa....)]}>AAA.A",
  ".....((((((.....))))))....(((((((((((....))))))))))).......(((((((((((.[[[[[)))))))))))((((((((((((....))))).)))))))....]]]]].(((((((((((((....))))))))))))).....................",
  ".....((((((.....))))))....(((((((((((....)))))))))))........((((((((((.[[[[[))))))))))(((((((((((((....))))).))))))))...]]]]].(((((((((((((....))))))))))))).....................",
  "((.(.(.((.(.....).)).)....(((((((((((....)))))))))))......)((((((((((...{{{{.)))))))))).(((((((((((....))))).))))))...)[}}}}..[[[[[[[[[[[[[)...]]]]]]]]]]]]]...].................",
  "(((((((([.{)](..[[.{)].](.}).})))))))",
  ".(((((..(((((....)))))((((((((.((........((((((..(((((((..)))))))(((((([..{{{{{{{)))))))))))..).))).))))))))))))].}}}}}}}...",
  "((((([{<aa<....))))).......]}>>AA",
  ".....(.((.(.....).)).)....(.((..({(<.{.a.a{........(..(((((.(.[....).)))).)...)).).}).((.(}...}...)>)).A.A...).)..............(((((((((((.(....).)))))))))))............]........",
  "..((...((.........))......(((((((....))))))).(((((((....))))))).(((((((....))))))).(((((((((....)))))))))..(((((((....)))))))..(((((((....)))))))...[[[[[[[......[[[[[[[[......[[.....)).....((..........]]......]]]]]]]]......................))..............]]]]]]].....(((((((....)))))))......................",
  "((((((((((((((((....)))))))((((.(((((((((((......)))))[[[[[[.)))))).))))(((((((((((((.((.(.....(......(((.(((((((((((.((((.(((((..]]]]]].))))).)))).))(((((((....)))))))(((((((((....)))))))(((((((....)))))))(((.(.((((((((((......)))))[[[[[[.))))).).)))))))))))))).)))..............)..).)).)))))))))))(((((((.(((((..]]]]]].))))).)))))))(((((((....))))))))))))))))))",
  "(((....))).(((....)))",
  "(((....))).(((....)))",
]

def test_structure_sanitization_success():

  for (i, case) in enumerate(success_cases):
    bp_list = utils.convert_dotbracket_to_bp_list(case[0], allow_pseudoknots=case[1])
    dbn = utils.convert_bp_list_to_dotbracket(bp_list, seq_len=len(case[0]))
    assert(dbn == success_expected_output[i])

failure_cases = [
  ["(((...))))", False],
  ["(((...))))", True],
  ["(((", False],
  ["(((", True],
  ["...)))", False],
  ["...)))", True],
  ["xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx", False],
  ["xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx", True],
  ["(((.[[[.(((...))).]]].)))", False],
  ["aaa.....AAA", False]
]
failure_expected_output = [
  "Unbalanced parenthesis notation: found closing character ')'",
  "Unbalanced parenthesis notation: found closing character ')'",
  "Unbalanced parenthesis notation: found unclosed pair for character '('",
  "Unbalanced parenthesis notation: found unclosed pair for character '('",
  "Unbalanced parenthesis notation: found closing character ')'",
  "Unbalanced parenthesis notation: found closing character ')'",
  "Unexpected character 'x'; did you mean to pass allow_pseudoknots=True?",
  "Unbalanced parenthesis notation: found unclosed pair for character 'x'",
  "Mixed pair delimiters found: '[' and '('; did you mean to pass allow_pseudoknots=True?",
  "Unexpected character 'a'; did you mean to pass allow_pseudoknots=True?"
]

def test_structure_sanitization_failure():
  for (i, case) in enumerate(failure_cases):
    with pytest.raises(Exception) as exc_info:
      bp_list = utils.convert_dotbracket_to_bp_list(case[0], allow_pseudoknots=case[1])
    assert(str(exc_info.value) == failure_expected_output[i])