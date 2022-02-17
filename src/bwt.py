"""Implementatin of the Burrows-Wheeler transform and related algorithms."""

from typing import (
    Iterator, Callable,
    NamedTuple
)

from alphabet import Alphabet
from sais import sais_alphabet
from subseq import SubSeq

ExactSearchFunc = Callable[[str], Iterator[int]]


def burrows_wheeler_transform(
    x: str
) -> tuple[bytearray, Alphabet, list[int]]:
    """
    Construct the Burrows-Wheeler transform.

    Build the bwt string from a string x, first mapping
    x to a new alphabet.

    Returns the transformed string, the mapping alphabet,
    and the suffix array over x.
    """
    x_, alpha = Alphabet.mapped_string_with_sentinel(x)
    sa = sais_alphabet(SubSeq[int](x_), alpha)
    bwt = bytearray(x_[j - 1] for j in sa)
    return bwt, alpha, sa


class CTable:
    """
    C-table for other bwt/fm-index search algorithms.

    for CTable ctab, ctab[⍺] is the number of occurrences
    of letters a < ⍺ in the bwt string (or the orignal string,
    since they have the same letters).
    """

    _cumsum: list[int]

    def __init__(self, bwt: bytearray, asize: int) -> None:
        """
        Construct a C-table.

        Compute the C-table from the bwt transformed string and
        the alphabet size.
        """
        # Count occurrences of characters in bwt
        counts = [0] * asize
        for a in bwt:
            counts[a] += 1
        # Get the cumulative sum
        n = 0
        for a, count in enumerate(counts):
            counts[a] = n
            n += count
        # That is all we need...
        self._cumsum = counts

    def __getitem__(self, a: int) -> int:
        """Get the number of occurrences of letters in the bwt less than a."""
        return self._cumsum[a]


class OTable:
    """
    O-table for the FM-index based search.

    For OTable otab, otab[a,i] is the number of occurrences j < i
    where bwt[j] == a.
    """

    _tbl: list[list[int]]

    def __init__(self, bwt: bytearray, asize: int) -> None:
        """
        Create O-table.

        Compute the O-table from the bwt transformed string and the size
        of the alphabet the bwt string is over.
        """
        # We exclude $ from lookups, so there are this many
        # rows.
        nrow = asize - 1
        # We need to index to len(bwt), but we don't represent first column
        # so there are len(bwt) columns.
        ncol = len(bwt)

        self._tbl = [[0] * ncol for _ in range(nrow)]

        # The first column is all zeros, the second
        # should hold a 1 in the row that has character
        # bwt[0]. The we b-1 because of the sentinel and
        # we use column 0 for the first real column.
        self._tbl[bwt[0] - 1][0] = 1

        # We already have cols 0 and 1. Now we need to
        # go up to (and including) len(bwt).
        for i in range(2, len(bwt) + 1):
            b = bwt[i - 1]
            # Characters, except for sentinel
            for a in range(1, asize):
                self._tbl[a - 1][i - 1] = self._tbl[a - 1][i - 2] + (a == b)

    def __getitem__(self, idx: tuple[int, int]) -> int:
        """
        Get the number of occurrences j < i where bwt[j] == a.

        a is the first and i the second value in the idx tuple.
        """
        a, i = idx
        assert a > 0, "Don't look up the sentinel"
        return 0 if i == 0 else self._tbl[a - 1][i - 1]


class FMIndexTables(NamedTuple):
    """Preprocessed FMIndex tables."""

    alpha: Alphabet
    sa: list[int]
    ctab: CTable
    otab: OTable


def preprocess_exact(x: str) -> FMIndexTables:
    """Preprocess tables for exact FM/bwt search."""
    bwt, alpha, sa = burrows_wheeler_transform(x)
    ctab = CTable(bwt, len(alpha))
    otab = OTable(bwt, len(alpha))
    return FMIndexTables(alpha, sa, ctab, otab)


def exact_searcher_from_tables(tbls: FMIndexTables) -> ExactSearchFunc:
    """Build an exact search function from preprocessed tables."""
    alpha, sa, ctab, otab = tbls

    def search(p_: str) -> Iterator[int]:
        try:
            p = alpha.map(p_)
        except KeyError:
            return  # can't map, so no matches

        # Find interval of matches...
        left, right = 0, len(sa)
        for a in reversed(p):
            left = ctab[a] + otab[a, left]
            right = ctab[a] + otab[a, right]
            if left >= right:
                return  # no matches

        # Report the matches
        for i in range(left, right):
            yield sa[i]

    return search


def exact_preprocess(x: str) -> ExactSearchFunc:
    """Build an exact search function for searching in string x."""
    return exact_searcher_from_tables(preprocess_exact(x))
