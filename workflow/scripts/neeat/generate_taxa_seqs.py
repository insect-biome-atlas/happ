import subprocess
import pandas as pd
from argparse import ArgumentParser
import os


def main(args):
    os.makedirs(args.outdir, exist_ok=True)
    # read in the taxonomy file
    taxdf = pd.read_csv(args.taxfile, sep="\t", index_col=0)
    # extract representative ASVs
    taxdf = taxdf.loc[taxdf["representative"] == 1]
    # also make sure that only clusters that remain in the counts file are used
    clusters_w_counts = []
    with open(args.countsfile, "r") as fhin:
        for line in fhin:
            clusters_w_counts.append(line.split("\t")[0])
    taxdf = taxdf.loc[taxdf["cluster"].isin(clusters_w_counts)]
    # taxa is the unique set of values for split_rank
    taxa = taxdf[args.rank].unique()
    singles = pd.DataFrame()
    # iterate over taxa
    for tax in taxa:
        # get asvs and cluster designation
        asvs = taxdf.loc[taxdf[args.rank] == tax, "cluster"]
        # skip taxa with fewer than 2 ASVs
        if len(asvs) < 2:
            singles = pd.concat([singles, asvs])
            continue
        # generate strings with '<asv> <cluster>' format to match fasta headers
        asvs = [f"{x[0]} {x[1]}" for x in list(zip(asvs.index, asvs))]
        # write asvs to file
        with open(f"{args.outdir}/{tax}.ids", "w") as f:
            _ = f.write("\n".join(asvs))
        # use seqkit to extract sequences
        with open(f"{args.outdir}/{tax}.fasta", "w") as fhout:
            cmd = [
                "seqkit",
                "grep",
                "-f",
                f"{args.outdir}/{tax}.ids",
                "-n",
                "--quiet",
                args.fastafile,
            ]
            _ = subprocess.run(cmd, stdout=fhout)
            os.remove(f"{args.outdir}/{tax}.ids")
    singles.to_csv(f"{args.outdir}/single_otus.tsv", sep="\t")


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
