"""Code for preprocessing a genome for FM-index search."""

import pickle
from bwt import (
    preprocess_exact,
    exact_searcher_from_tables,
    ExactSearchFunc
)


# A genome maps from chromosome names to chromosome sequences
GENOME = dict[str, str]
# Once preprocessed and loaded, we have a table of search functions
# instead.
GENOME_SEARCH = dict[str, ExactSearchFunc]


def preprocess(genome: GENOME, preproc_file_name: str) -> None:
    """Preprocess a genome and pickle the generated search functions."""
    preprocessed = {
        name: preprocess_exact(seq) for name, seq in genome.items()
    }
    with open(preproc_file_name, "wb") as preproc_file:
        pickle.dump(preprocessed, preproc_file)


def load_preprocessed(preproc_file_name: str) -> GENOME_SEARCH:
    """Load preprocessed tables and make them into search functions."""
    with open(preproc_file_name, "rb") as preproc_file:
        preproc_tables = pickle.load(preproc_file)
    return {
        name: exact_searcher_from_tables(tbl)
        for name, tbl in preproc_tables.items()
    }
