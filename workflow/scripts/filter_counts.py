#!/usr/bin/env python

from argparse import ArgumentParser
import polars as pl
from Bio.SeqIO import parse
import logging
import gzip as gz


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(message)s",
)


def read_records(fastafile):
    records = []
    if fastafile.endswith(".gz"):
        open_fn = gz.open
    else:
        open_fn = open
    with open_fn(fastafile, "rt") as fhin:
        for record in parse(fhin, "fasta"):
            records.append(record.id)
    return records


def main(args):
    counts = pl.scan_csv(args.counts, separator="\t")
    idcol = counts.collect_schema().names()[0]
    records = read_records(args.fasta)
    outfile = f"{args.outdir}/asv_counts.tsv"
    counts.filter(pl.col(idcol).is_in(records)).sink_csv(outfile, separator="\t")
    total_counts = (
        counts.filter(pl.col(idcol).is_in(records))
        .with_columns(total=pl.sum_horizontal(pl.col("*").exclude(idcol)))
        .select(idcol, "total")
        .collect()
    )
    total_counts.columns = ["Representative_Sequence", "total"]
    total_counts = total_counts.select(
        pl.col("Representative_Sequence"), pl.col("total").cast(pl.Int32)
    )
    outfile = f"{args.outdir}/total_counts.tsv"
    total_counts.write_csv(outfile, separator="\t")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-f", "--fasta", help="Input fasta file")
    parser.add_argument("-c", "--counts", help="Input counts file")
    parser.add_argument("-o", "--outdir", help="Output directory")
    args = parser.parse_args()
    main(args)
