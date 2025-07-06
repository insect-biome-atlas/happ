import subprocess
import pandas as pd
from argparse import ArgumentParser
import os
import sys

def check_fasta_header(f):
    with open(f, 'r') as fhin:
        line = fhin.readline()
    return line.rstrip().lstrip(">").split(" ")

def main(args):
    os.makedirs(args.outdir, exist_ok=True)
    # read in the taxonomy file
    taxdf = pd.read_csv(args.taxfile, sep="\t", index_col=0)
    # extract representative sequences
    if "representative" in taxdf.columns:
        taxdf = taxdf.loc[taxdf.representative==1]
    taxdf.index.name="seqid"
    if "cluster" in taxdf.columns:
        taxdf = taxdf.reset_index().set_index("cluster")
    # also make sure that only clusters that remain in the counts file are used
    clusters_w_counts = []
    with open(args.countsfile, "r") as fhin:
        clusters_w_counts = [line.split("\t")[0] for line in fhin]
    taxdf = taxdf.loc[taxdf.index.isin(clusters_w_counts)]
    # taxa is the unique set of values for split_rank
    taxa = taxdf[args.rank].unique()
    singles = pd.DataFrame()
    header = check_fasta_header(args.fastafile)
    # if fasta header only contains one identifier
    # we assume this identifier is also in the taxdf index
    if len(header) == 1 and header[0] in taxdf.index:
        cluster_desc = False
    # otherwise we assume the second identifier in the header matches
    # the taxdf index
    elif len(header) > 1 and header[1] in taxdf.index:
        cluster_desc = True
    else:
        sys.exit("Cannot match fasta headers to taxfile")
    # iterate over taxa
    for tax in taxa:
        # get asvs and cluster designation
        seqs = taxdf.loc[taxdf[args.rank] == tax]
        # skip taxa with fewer than 2 ASVs
        if len(seqs) < 2:
            singles = pd.concat([singles, seqs])
            continue
        cmd = f"seqkit grep -f {args.outdir}/{tax}.ids -n --quiet {args.fastafile}"
        seqs = seqs.index.tolist()
        if cluster_desc:
            # generate strings with '<asv> <cluster>' format to match fasta headers
            seqs = [f"{x[0]} {x[1]}" for x in list(zip(taxdf.loc[seqs, "seqid"], seqs))]
            # add seqkit replace command to strip first id from fasta headers
            cmd += "| seqkit replace -p '^.+ '"
        # write asvs to file
        with open(f"{args.outdir}/{tax}.ids", "w") as f:
            _ = f.write("\n".join(seqs))
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        output, _ = process.communicate()
        # use seqkit to extract sequences
        with open(f"{args.outdir}/{tax}.fasta", "w") as fhout:
            fhout.write(output.decode())
        os.remove(f"{args.outdir}/{tax}.ids")
    if len(singles) > 0:
        singles.to_csv(f"{args.outdir}/singles.tsv", sep="\t")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-t", "--taxfile", type=str, help="Taxonomy file")
    parser.add_argument(
        "-f", "--fastafile", type=str, help="Fasta file of cluster/OTU representatives"
    )
    parser.add_argument(
        "-c",
        "--countsfile",
        type=str,
        help="TSV file with counts for clusters/OTUs (rows) in samples (columns)",
    )
    parser.add_argument(
        "-r", "--rank", type=str, help="Rank at which to partition sequences"
    )
    parser.add_argument(
        "-o", "--outdir", type=str, help="Write sequences to files in outdir"
    )
    args = parser.parse_args()
    main(args)
