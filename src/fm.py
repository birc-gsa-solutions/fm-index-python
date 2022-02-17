"""FM-index exact pattern matching."""

import argparse
import sys

from preprocess import (
    preprocess,
    load_preprocessed
)
from fasta import read_fasta
from fastq import scan_reads
from sam import ssam_record


def main() -> None:
    """FM-index exact pattern matching."""
    argparser = argparse.ArgumentParser(
        description="FM-index exact pattern matching",
        usage="\n\tfm -p genome\n\tfm genome reads"
    )
    argparser.add_argument(
        "-p", action="store_true",
        help="preprocess the genome."
    )
    argparser.add_argument(
        "genome",
        help="Simple-FASTA file containing the genome.",
        type=argparse.FileType('r')
    )
    argparser.add_argument(
        "reads", nargs="?",
        help="Simple-FASTQ file containing the reads.",
        type=argparse.FileType('r')
    )
    args = argparser.parse_args()

    if args.p:
        preprocess(read_fasta(args.genome), args.genome.name+".fm-idx")
    else:
        # here we need the optional argument reads
        if args.reads is None:
            argparser.print_help()
            sys.exit(1)

        genome_searchers = load_preprocessed(args.genome.name+".fm-idx")
        for read_name, read_seq in scan_reads(args.reads):
            for chr_name, search in genome_searchers.items():
                for i in search(read_seq):
                    ssam_record(sys.stdout,
                                read_name, chr_name,
                                i, f"{len(read_seq)}M",
                                read_seq)


if __name__ == '__main__':
    main()
