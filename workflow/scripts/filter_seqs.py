#!/usr/bin/env python

from argparse import ArgumentParser
import polars as pl
from Bio.SeqIO import parse
import subprocess
from tqdm import tqdm
import os
import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(message)s",
)


def rename_records(fastafile, rename_dict, outfile):
    records = []
    with open(fastafile, "r") as fhin, open(outfile, "w") as fhout:
        for record in parse(fhin, "fasta"):
            try:
                newid = rename_dict[record.id]
            except KeyError:
                continue
            records.append(record.id)
            fhout.write(f">{newid}\n{record.seq}\n")
    return records


def main(args):
    rank = args.rank
    logging.info(f"Reading taxonomy from {args.taxfile}")
    tax = pl.read_csv(args.taxfile, separator="\t")
    taxa = tax.select(rank).unique().to_series().to_list()
    logging.info(f"{len(taxa)} unique taxa at {args.rank}")
    idcol = tax.columns[0]
    tax = tax.select(idcol, rank)
    asvs = tax.select(idcol).to_series().to_list()
    renames = (
        tax.with_columns(rename=pl.concat_str(idcol, rank, separator=";"))
        .select("rename")
        .to_series()
        .to_list()
    )
    rename_dict = dict(zip(asvs, renames))
    os.makedirs(args.tmpdir, exist_ok=True)
    tempfile = f"{args.tmpdir}/renamed.fasta"
    records = rename_records(args.fastafile, rename_dict, tempfile)
    logging.info(f"{len(records)} unique ASVs found in {args.fastafile}")
    splitdir = f"{args.tmpdir}/splits"
    cmd = [
        "seqkit",
        "split",
        "--force",
        "-j",
        str(args.threads),
        "-i",
        "--id-regexp",
        ";(.+)",
        "--quiet",
        "-O",
        splitdir,
        tempfile,
    ]
    returncode = subprocess.run(cmd)
    for t in tqdm(taxa, desc=f"partitioning taxa by {rank}", unit=f" taxa", ncols=120):
        outdir = f"{args.outdir}/{t}"
        t_fasta = f"{splitdir}/renamed.part_{t}.fasta"
        if not os.path.exists(t_fasta):
            continue
        os.makedirs(outdir, exist_ok=True)
        t_fasta_out = f"{outdir}/asv_seqs.fasta.gz"
        cmd = [
            "seqkit",
            "replace",
            "-j",
            str(args.threads),
            "-p",
            "(.+);.+",
            "-r",
            "$1",
            "-o",
            t_fasta_out,
            t_fasta,
        ]
        returncode = subprocess.run(cmd)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-t", "--taxfile", help="Taxonomy assignments file")
    parser.add_argument("-f", "--fastafile", help="Fasta file")
    parser.add_argument("-r", "--rank", default="Family", help="Rank to split by")
    parser.add_argument("-o", "--outdir", help="Output directory", default="out")
    parser.add_argument("-j", "--threads", help="Number of CPUs", default=1)
    parser.add_argument("--tmpdir", help="Temporary directory", default="temp")
    args = parser.parse_args()
    main(args)
