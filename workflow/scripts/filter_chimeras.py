#!/usr/bin/env python

from argparse import ArgumentParser
import pandas as pd
from Bio.SeqIO import parse
import sys
import gzip as gz
import numpy as np
import time



def add_size_column(df):
    size = [int(i.split("=")[-1]) for i in df.index]
    df.rename(index=lambda x: x.split(";")[0], inplace=True)
    return df.assign(Size=pd.Series(size, index=df.index))


def filter_uchime(df, minh, mindiff, mindiv):
    """
    Performs first round filtering according to uchime settings
    """
    return df.loc[(df.Score>=minh)&(df.sumL>=mindiff)&(df.sumR>=mindiff)&(df.Div>=mindiv)]


def read_uchime(f):
    """
    Columns:
    0. Score: Chimera score
    1. Q: Query label.
    2. A: Parent A label.
    3. B: Parent B label.
    4. T: Top parent (T) label. This is the closest reference sequence; usually either A or B.
    5. IdQM (QM): Percent identity of query and the model (M) constructed as a segment of A and a segment of B.
    6. IdQA (QA): Percent identity of Q and A.
    7. IdQB (QB): Percent identity of Q and B.
    8. IdAB (AB): Percent identity of A and B
    9. IdQT (QT): Percent identity of Q and T.
    10. LY (best_left_y): Yes votes in left segment.
    11. LN (best_left_n): No votes in left segment.
    12. LA (best_left_a): Abstain votes in left segment.
    13. RY (best_right_y): Yes votes in right segment.
    14. RN (best_right_n): No votes in right segment.
    15. RA (best_right_a): Abstain votes in right segment.
    16. Div (divdiff): Divergence, defined as (IdQM - IdQT).
    17. YN (status): Y or N, indicating whether the query was classified as chimeric. This requires that Score >= threshold specified by -minh, Div > minimum divergence specified by ‑mindiv and the number of diffs ( (Y+N+A) in each segment (L and R) is greater than the minimum specified by -mindiffs. In v6.0.310 and later, may also be '?' indicating a weakly chimeric alignment with score between maxh and minh.
    """
    names = [
        "Score",
        "Q",
        "A",
        "B",
        "T",
        "IdQM",
        "IdQA",
        "IdQB",
        "IdAB",
        "IdQT",
        "LY",
        "LN",
        "LA",
        "RY",
        "RN",
        "RA",
        "Div",
        "YN",
    ]
    df = pd.read_csv(f, sep="\t", index_col=1, names=names)
    df = add_size_column(df)
    df = df.replace("*", np.nan)
    for col in ["IdQM", "IdQA", "IdQB", "IdAB", "IdQT", "Div", "Size"]:
        df[col] = pd.to_numeric(df[col])
    df["sumL"] = df["LN"] + df["LA"] + df["LY"]
    df["sumR"] = df["RN"] + df["RA"] + df["RY"]
    return df


def calc_occurrence(f, chunksize=1000):
    start = time.perf_counter()
    d = {}
    n_samples = 0
    with pd.read_csv(f, sep="\t", index_col=0, chunksize=chunksize) as reader:
        reader
        for i, chunk in enumerate(tqdm(reader, unit=" chunks",
                          desc="calculating occurrence")):
            d.update(chunk.apply(lambda x: x.index[x > 0].to_list(),
                                 axis=1).to_dict())
            if i==0:
                n_samples = chunk.shape[1]
    stop = time.perf_counter()
    n_asvs = len(list(d.keys()))
    sys.stderr.write(f"Read {n_asvs} ASVs across {n_samples} samples in {stop - start} seconds\n")
    return d, n_samples


def add_shared_cols(df, d, col, n_samples):
    n_shared_vals = []
    frac_shared_vals = []
    for i, q in enumerate(
            tqdm(df.index, unit=" ASVs", desc=f"getting shared samples for {col}")
    ):
        p = df[col][i].split(";")[0]
        n_shared, frac_shared = shared(d, q, p, n_samples)
        n_shared_vals.append(n_shared)
        frac_shared_vals.append(frac_shared)
    df = df.assign(shared=pd.Series(n_shared_vals, index=df.index))
    df = df.assign(frac_shared=pd.Series(frac_shared_vals, index=df.index))
    return df.rename(columns={'shared': f"{col}_shared", 'frac_shared': f"{col}_frac_shared"})


def shared(d, asv1, asv2, n_samples):
    asv1_samples = len(set(d[asv1]))
    n_shared = len(set(d[asv1]).intersection(d[asv2]))
    frac_shared = n_shared / asv1_samples
    return n_shared, frac_shared


def main(args):
    if len(args.uchimeout) > 1:
        # Run samplewise chimera detection
    else:
        sys.stderr.write("Only one chimera results file found. Running in batchwise mode\n")
        sys.stderr.write(f"Reading chimera results from {args.uchimeout[0]}\n")
        uchime_res = read_uchime(args.uchimeout)
        uchime_filtered = filter_uchime(uchime_res, args.minh, args.mindiff, args.mindiv)
        if args.min_samples_shared > 0 or args.min_frac_samples_shared > 0:
            min_samples_shared = args.min_samples_shared
            min_frac_samples_shared = args.min_frac_samples_shared
            d, n_samples = calc_occurrence(args.counts, chunksize=args.chunksize)
            # Add number of shared samples + shared fraction for each Q->A, and Q->B pair
            uchime_filtered = add_shared_cols(uchime_filtered, d, "A", n_samples)
            uchime_filtered = add_shared_cols(uchime_filtered, d, "B", n_samples)
            filtered = uchime_filtered.loc[(uchime_filtered.A_shared>=min_samples_shared)&(uchime_filtered.A_frac_shared>=min_frac_samples_shared)&(uchime_filtered.B_shared>=min_samples_shared)&(uchime_filtered.B_frac_shared>=min_frac_samples_shared)]
            chimeras = list(filtered.index)
        else:
            chimeras = list(uchime_filtered.index)



if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--counts", type=str, help="ASV counts file")
    parser.add_argument("--fasta", type=str, help="ASV fasta file")
    parser.add_argument("--uchimeout", type=str, nargs="+", help="Uchime results file(s)")
    parser.add_argument("--min_samples_shared", type=int, default=0,
                        help="Minimum number of samples in which a sequence has "
                             "to be present with its parents for it to be classified"
                             "as a chimera")
    parser.add_argument("--min_frac_samples_shared", type=float, default=0,
                        help="Minimum fraction of samples in which a sequence "
                             "has to be present with its parents for it to be"
                             "classified as a chimera")
    parser.add_argument("--min_chimeric_samples", type=int,
                        help="Minimum number of samples in which a sequence"
                             "has to be marked as chimeric for it to be classified as a chimera")
    parser.add_argument("--mindiff", type=float, default=3,
                        help="Minimum number of diffs in a segment.")
    parser.add_argument("--mindiv", type=float, default=0.8,
                        help="Minimum divergence, i.e. 100% - identity between "
                             "the query and closest reference database sequence.")
    parser.add_argument("--minh", type=float, default=0.28,
                        help="Minimum score (h) to be classified as chimera.")
    parser.add_argument("--chunksize", type=int, default=1000,
                        help="Read counts file in chunks of <chunk> rows")
    args = parser.parse_args()
    main(args)
