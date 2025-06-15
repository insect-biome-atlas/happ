#!/usr/bin/env python

from argparse import ArgumentParser, SUPPRESS
import pandas as pd
from Bio.SeqIO import parse
import tqdm
import logging
import sys
from tqdm.contrib.logging import logging_redirect_tqdm

def read_seqs(files, outfile):
    seqids = {}
    with open(outfile, 'w') as fhout:
        for file in files:
            logging.info(f"Reading sequences from {file}")
            start_seqs = len(seqids)
            for record in parse(file, "fasta"):
                if not record.id in seqids:
                    fhout.write(f">{record.id}\n{record.seq}\n")
                    seqids[record.id] = True
            end_seqs = len(seqids)
            logging.info(f"Wrote {end_seqs - start_seqs} sequences ({len(seqids)} total)")
    return seqids

def read_tax_tables(files, outfile, seqids):
    df = pd.DataFrame()
    for file in files:
        start = df.shape[0]
        logging.info(f"Reading tax table from {file}")
        _df = pd.read_csv(file, sep="\t", index_col=0)
        df = pd.concat([
            _df.loc[_df.index.isin(seqids)&(~_df.index.isin(df.index))], 
            df])
        end = df.shape[0]
        logging.info(f"Added {end - start} sequences to tax table ({df.shape[0]} total)")
    df.to_csv(outfile, sep="\t")
    

def read_count_tables(files, outfile, seqids, trim_sample_names=False, nrows=None):
    d = {}   
    for file in files:
        logging.info(f"Reading count table from {file}")
        with open(file, 'r') as fhin:
            pbar = tqdm.tqdm(fhin, unit=" lines", leave=True)
            i = -1
            for line in pbar:
                i+=1
                pbar.set_postfix({'dict_mem_size (GB)': sys.getsizeof(d)})
                items = line.rstrip().split("\t")
                if i == 0:
                    samples = items[1:]
                    if trim_sample_names:
                        samples = ["_".join(s.split("_")[0:2]) for s in samples]
                    continue
                asv = items[0]
                if asv not in seqids:
                    continue
                for j, item in enumerate(items[1:]):
                    count = int(float(item))
                    sample = samples[j]
                    if sample in d and asv in d[sample] and d[sample][asv] >= count:
                            continue
                    if sample not in d:
                        d[sample] = {}
                    d[sample][asv] = count
                if nrows and i >= nrows:
                    break
    counts = pd.DataFrame(d).fillna(0).astype(int)
    counts.to_csv(outfile, sep="\t")
        

def main(args):
    asv_count_tables = args.asv_count_tables
    asv_fasta_files = args.asv_fasta_files
    asv_tax_tables = args.asv_tax_tables
    outdir = args.outdir
    logging.basicConfig(filename=f"{outdir}/log", filemode="w", format='%(levelname)s:%(message)s', encoding='utf-8', level=logging.DEBUG)
    if not asv_fasta_files:
        logging.error("No ASV fasta files provided")
        sys.exit(1)
    # Read ASV fasta files
    seqids = read_seqs(asv_fasta_files, f"{outdir}/asv_seqs.fasta")
    
    if asv_tax_tables:
        # Read ASV tax tables
        read_tax_tables(asv_tax_tables, f"{outdir}/asv_taxa.tsv", seqids)

    if asv_count_tables:
        # Read ASV count tables
        read_count_tables(asv_count_tables, f"{outdir}/asv_counts.tsv", seqids, args.trim_sample_names, nrows=args.nrows)

if __name__ == "__main__":
    parser = ArgumentParser(description="Combine results from different tools")
    parser.add_argument("-c", "--asv_count_tables", nargs="+", help="ASV count tables")
    parser.add_argument("-f", "--asv_fasta_files", nargs="+", help="ASV fasta files")
    parser.add_argument("-t", "--asv_tax_tables", nargs="+", help="ASV taxonomy tables")
    parser.add_argument("-o", "--outdir", type=str, help="Output directory", required=True)
    parser.add_argument("--trim_sample_names", action="store_true", help="Trim sample names to first two parts separated by underscores")
    parser.add_argument("--nrows", type=int, help=SUPPRESS, default=None)
    args = parser.parse_args()
    main(args)