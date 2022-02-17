"""Test bwt."""

from test_helpers import check_matches
import alphabet
import bwt
import sais


def test_ctable() -> None:
    """Test C-table."""
    x, alpha = alphabet.Alphabet.mapped_string_with_sentinel("aabca")
    ctab = bwt.CTable(x, len(alpha))
    assert ctab[0] == 0, "Nothing is smaller than the sentinel"
    assert ctab[1] == 1, "$ is smaller than 'a'"
    assert ctab[2] == 4, "$ + three 'a'"
    assert ctab[3] == 5, "$ + three 'a' + one 'b'"


def test_otable() -> None:
    """Test O-table."""
    x = "aabca"
    transformed, alpha, _ = bwt.burrows_wheeler_transform(x)
    assert transformed == bytearray([1, 3, 0, 1, 1, 2])

    # we shouldn't look at private members, of course, but
    # we are only testing...
    # pylint: disable=protected-access
    otab = bwt.OTable(transformed, len(alpha))
    assert len(otab._tbl) == len(alpha) - 1
    assert len(otab._tbl[0]) == len(transformed)
    assert otab._tbl[0] == [1, 1, 1, 2, 3, 3], "a counts"
    assert otab._tbl[1] == [0, 0, 0, 0, 0, 1], "b counts"
    assert otab._tbl[2] == [0, 1, 1, 1, 1, 1], "c counts"


def test_mississippi() -> None:
    """Test on mississippi."""
    x = "mississippi"
    sa = sais.sais(x)
    for j in sa:
        print(x[j:])
    search = bwt.exact_preprocess(x)
    for p in ("si", "ppi", "ssi", "pip", "x", ""):
        matches = list(search(p))
        print(p, matches)
        check_matches(x, p, matches)
    # the empty string should give us the entire x includng sentinel
    assert len(list(search(""))) == len(x) + 1
