#!/usr/bin/env python

from argparse import ArgumentParser
from collections import defaultdict
import pandas as pd
from Bio.SeqIO import parse
import sys
from tqdm import tqdm
import numpy as np
import time


def add_size_column(df):
    size = [int(float(i.split("=")[-1])) for i in df.index]
    df.rename(index=lambda x: x.split(";")[0], inplace=True)
    return df.assign(Size=pd.Series(size, index=df.index))


def filter_uchime(df, minh, mindiff, mindiv, algorithm):
    """
    Performs first round filtering according to uchime settings
    """
    if algorithm in ["uchime2_denovo", "uchime3_denovo"]:
        return df.loc[
            (df.IdQM == 100)
            & (df.IdQT == 100)
        ]
    else:
        return df.loc[
            (df.Score >= minh)
            & (df.sumL >= mindiff)
            & (df.sumR >= mindiff)
            & (df.Div >= mindiv)
        ]


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
    if df.shape[0] == 0:
        return df
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
        for i, chunk in enumerate(
            tqdm(reader, unit=" chunks", desc="calculating occurrence")
        ):
            d.update(chunk.apply(lambda x: x.index[x > 0].to_list(), axis=1).to_dict())
            if i == 0:
                n_samples = chunk.shape[1]
    stop = time.perf_counter()
    n_asvs = len(list(d.keys()))
    sys.stderr.write(
        f"Read {n_asvs} ASVs across {n_samples} samples in {stop - start} seconds\n"
    )
    return d, n_samples


def add_shared_cols(df, d, col):
    n_shared_vals = []
    frac_shared_vals = []
    for i, q in enumerate(
        tqdm(df.index, unit=" ASVs", desc=f"getting shared samples for {col}")
    ):
        try:
            p = df[col][i].split(";")[0]
            n_shared, frac_shared = shared(d, q, p)
            n_shared_vals.append(n_shared)
            frac_shared_vals.append(frac_shared)
        except AttributeError:
            n_shared_vals.append(np.nan)
            frac_shared_vals.append(np.nan)
    df = df.assign(shared=pd.Series(n_shared_vals, index=df.index))
    df = df.assign(frac_shared=pd.Series(frac_shared_vals, index=df.index))
    return df.rename(
        columns={"shared": f"{col}_shared", "frac_shared": f"{col}_frac_shared"}
    )


def shared(d, asv1, asv2):
    asv1_samples = len(set(d[asv1]))
    n_shared = len(set(d[asv1]).intersection(d[asv2]))
    frac_shared = n_shared / asv1_samples
    return n_shared, frac_shared


def write_seqs(f, d):
    with open(f, "w") as fhout:
        for record, text in tqdm(d.items(), unit=" seqs", desc=f"writing to {f}"):
            fhout.write(text)


def filter_samplewise_chimeras(files, minh, mindiff, mindiv, min_chimeric_samples, min_frac_chimeric_samples, algorithm):
    all_seqs = []
    asvs = defaultdict(lambda: 0)
    asvs_sample_counts = defaultdict(lambda: 0)
    nonchims = []
    for uchimeout in tqdm(files, unit=" files", desc="Reading chimera results"):
        uchime_res = read_uchime(uchimeout)
        # If empty input, skip to next file
        if uchime_res.shape[0] == 0:
            continue
        all_asvs = set(uchime_res.index)
        # Count number of samples in which ASVs were found
        for asv in all_asvs:
            asvs_sample_counts[asv] += 1
        all_seqs = list(set(all_seqs + list(uchime_res.index)))
        uchime_filtered = filter_uchime(uchime_res, minh, mindiff, mindiv, algorithm)
        chim_asvs = set(uchime_filtered.index)
        asv_diffs = all_asvs.difference(chim_asvs)
        # if there are asvs in all_asvs that are not in chim_asvs and
        # min_chimeric_samples == 0 or min_frac_chimeric_samples == 0, 
        # set these asvs to nonchims
        # this means that as soon as an ASV is not identified as chimeric in a sample
        # it should not be treated as chimeric in other samples
        if min_chimeric_samples == 0 or min_frac_chimeric_samples == 0:
            nonchims += list(asv_diffs)
        # Keep track of number of times each ASV is marked as chimeric
        for asv in chim_asvs:
            asvs[asv] += 1
    chimcounts = pd.DataFrame(asvs, index=["n"]).T
    sample_counts = pd.DataFrame(asvs_sample_counts, index=["in_n_samples"]).T
    # Add column on fraction of samples in which an ASV was identified as chimeric
    chimcounts = pd.merge(chimcounts, sample_counts, left_index=True, right_index=True)
    chimcounts["frac"] = chimcounts["n"].div(chimcounts["in_n_samples"])
    if min_frac_chimeric_samples:
        chimeras = list(chimcounts.loc[chimcounts.frac >= min_frac_chimeric_samples].index)
    else:
        # Set chimeras to ASVs identified as chimeric in at least <min_chimeric_samples>
        chimeras = list(chimcounts.loc[chimcounts.n >= min_chimeric_samples].index)
    # However, remove ASVs already marked as non-chimeric
    chimeras = list(set(chimeras).difference(nonchims))
    nonchimeras = list(set(all_seqs).difference(chimeras))
    return chimeras, nonchimeras


def main(args):
    if len(args.uchimeout) > 1:
        # Run samplewise chimera detection
        sys.stderr.write(f"Found {len(args.uchimeout)} chimera result files\n")
        sys.stderr.write(
            f"Running in samplewise mode using: minh={args.minh}, mindiff={args.mindiff}, mindiv={args.mindiv}, min_chimeric_samples={args.min_chimeric_samples}, min_frac_chimeric_samples={args.min_frac_chimeric_samples}\n"
        )
        chimeras, nonchimeras = filter_samplewise_chimeras(
            args.uchimeout,
            args.minh,
            args.mindiff,
            args.mindiv,
            args.min_chimeric_samples,
            args.min_frac_chimeric_samples,
            args.algorithm,
        )
    else:
        sys.stderr.write(
            "Only one chimera results file found. Running in batchwise mode\n"
        )
        sys.stderr.write(f"Reading chimera results from {args.uchimeout[0]}\n")
        uchime_res = read_uchime(args.uchimeout[0])
        min_samples_shared = args.min_samples_shared
        min_frac_samples_shared = args.min_frac_samples_shared
        sys.stderr.write(f"Calculating occurrence based on {args.counts}\n")
        d, n_samples = calc_occurrence(args.counts, chunksize=args.chunksize)
        # Add number of shared samples + shared fraction for each Q->A, and Q->B pair
        uchime_res = add_shared_cols(uchime_res, d, "A")
        uchime_res = add_shared_cols(uchime_res, d, "B")
        # Filter according to uchime params
        sys.stderr.write(
            f"1st round of filtering with minh={args.minh}, mindiff={args.mindiff}, mindiv={args.mindiv}\n"
        )
        uchime_filtered = filter_uchime(
            uchime_res, args.minh, args.mindiff, args.mindiv, args.algorithm
        )
        sys.stderr.write(f"{uchime_filtered.shape[0]} potential chimeras remaining\n")
        # Filter according to shared samples
        if min_samples_shared > 0 or min_frac_samples_shared > 0:
            sys.stderr.write(
                f"2nd round of filtering with min_samples_share={min_samples_shared}, min_frac_samples_shared={min_frac_samples_shared}\n"
            )
            uchime_filtered = uchime_filtered.loc[
                (uchime_filtered.A_shared >= min_samples_shared)
                & (uchime_filtered.A_frac_shared >= min_frac_samples_shared)
                & (uchime_filtered.B_shared >= min_samples_shared)
                & (uchime_filtered.B_frac_shared >= min_frac_samples_shared)
            ]
            sys.stderr.write(
                f"{uchime_filtered.shape[0]} potential chimeras remaining\n"
            )
        chimeras = list(uchime_filtered.index)
        nonchimeras = set(uchime_res.index).difference(chimeras)
        if args.filteredout:
            uchime_res.to_csv(args.filteredout, sep="\t")
    sys.stderr.write(f"{len(nonchimeras)} non-chimeric seqs identified\n")
    if args.fasta:
        nonchimseqs = {}
        chimseqs = {}
        for record in tqdm(
            parse(args.fasta, "fasta"), unit=" seqs", desc="Reading ASV seqs"
        ):
            if args.chimfasta and record.id in chimeras:
                chimseqs[record.id] = f">{record.id}\n{record.seq}\n"
            elif args.nonchimfasta and record.id in nonchimeras:
                nonchimseqs[record.id] = f">{record.id}\n{record.seq}\n"
        if args.chimfasta:
            sys.stderr.write(
                f"Writing {len(chimseqs.keys())} chimeric seqs to {args.chimfasta}\n"
            )
            write_seqs(args.chimfasta, chimseqs)
        if args.nonchimfasta:
            sys.stderr.write(
                f"Writing {len(nonchimseqs.keys())} non-chimeric seqs to {args.nonchimfasta}\n"
            )
            write_seqs(args.nonchimfasta, nonchimseqs)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--counts", type=str, help="ASV counts file")
    parser.add_argument("--fasta", type=str, help="ASV fasta file")
    parser.add_argument("--chimfasta", type=str, help="Fasta file with chimeras")
    parser.add_argument("--nonchimfasta", type=str, help="Fasta file with nonchimeras")
    parser.add_argument(
        "--uchimeout",
        nargs="+",
        help="Uchime results file(s)",
        required=True,
    )
    parser.add_argument(
        "--filteredout", type=str, help="Write uchime results with additional info"
    )
    parser.add_argument(
        "--min_samples_shared",
        type=int,
        default=0,
        help="Minimum number of samples in which a sequence has to be present with its parents for it to be classified as a chimera",
    )
    parser.add_argument(
        "--min_frac_samples_shared",
        type=float,
        default=0.0,
        help="Minimum fraction of samples in which a sequence has to be present with its parents for it to be classified as a chimera",
    )
    parser.add_argument(
        "--min_chimeric_samples",
        type=int,
        default=1,
        help="Minimum number of samples in which a sequence has to be marked as chimeric for it to be classified as a chimera",
    )
    parser.add_argument(
        "--min_frac_chimeric_samples",
        type=float,
        help="Minimum fraction of samples in which a sequence has to be marked as chimeric for it to be classified as a chimera",
    )
    parser.add_argument(
        "--mindiff", type=int, default=3, help="Minimum number of diffs in a segment."
    )
    parser.add_argument(
        "--mindiv",
        type=float,
        default=0.8,
        help="Minimum divergence, i.e. 100%% - identity between the query and closest reference database sequence.",
    )
    parser.add_argument(
        "--minh",
        type=float,
        default=0.28,
        help="Minimum score (h) to be classified as chimera.",
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        default=1000,
        help="Read counts file in chunks of <chunk> rows",
    )
    parser.add_argument(
        "--algorithm",
        type=str,
        default="uchime_denovo",
        choices=["uchime_denovo", "uchime2_denovo", "uchime3_denovo"],
        help="Algorithm to filter by. Note than when 'uchime2_denovo' or 'uchime3_denovo' is specified this requires"
             "100%% matching between the alignment and the model (perfect chimeras). "
    )
    args = parser.parse_args()
    main(args)
