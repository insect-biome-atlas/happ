# /usr/bin/env python

from argparse import ArgumentParser
import pandas as pd
from Bio.SeqIO import parse
import sys
from collections import defaultdict
import numpy as np


def get_special_asvs(special_asv_seqs_file):
    asvid_isspecial = {}
    for record in parse(special_asv_seqs_file, "fasta"):
        asvid_isspecial[record.id] = 1
    return asvid_isspecial


def get_taxonomy(taxonomy_file):
    taxdf = pd.read_csv(taxonomy_file, sep="\t", header=0, index_col=0)
    return taxdf


def add_size_column(df):
    size = [int(i.split("=")[-1]) for i in df.index]
    df.rename(index=lambda x: x.split(";")[0], inplace=True)
    return df.assign(Size=pd.Series(size, index=df.index))


def count_otus(df):
    return len(df["cluster"].unique())


def count_species(df):
    return len(df["Species"].unique())


def read_uchime(f, nonchims=[]):
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
    df["Type"] = "NOF"
    df.loc[list(set(nonchims).intersection(df.index)), "Type"] = "FIN"
    df["sumL"] = df["LN"] + df["LA"] + df["LY"]
    df["sumR"] = df["RN"] + df["RA"] + df["RY"]
    return df


def eval_uchime(
    uchime_res,
    asv_taxonomy,
    minh="0.001,0.005,0.01,0.02,0.04,0.08,0.16,0.28,0.5,1.0,1000",
    mindiffs="0,1,2,3,4",
    mindiv="0.5,0.8,1.2,1.5",
):
    minh = [float(x) for x in minh.split(",")]
    mindiffs = [int(x) for x in mindiffs.split(",")]
    mindiv = [float(x) for x in mindiv.split(",")]
    uchime_res = pd.merge(uchime_res, asv_taxonomy, left_index=True, right_index=True)
    with sys.stdout as fhout:
        fhout.write(
            "minh\tASVs\totus\tbins\tspecies\treads\tFIN_ASVs\tFIN_otus\tFIN_bins\tFIN_species\tFIN_reads\tzero_bin_otus\tzero_species_otus\tmulti_otu_species\tzero_bin_otus_lepidoptera\tzero_species_otus_lepidoptera\tmulti_otu_species_lepidoptera\n"
        )
        for h in minh:
            for mindiff in mindiffs:
                for m in mindiv:
                    sys.stderr.write(f"#minh:{h}, mindiff:{mindiff}, mindiv:{m}\n")
                    res = uchime_res.loc[uchime_res.Score < h]
                    #((divdiff >= opt_mindiv) && (sumL >= opt_mindiffs) && (sumR >= opt_mindiffs))
                    res = res.loc[(res.Div<m)|(res.sumL<mindiff)|(res.sumR<mindiff)]
                    ## FINBOL results
                    fin_res = res.loc[res.Type == "FIN"]
                    ## Lepidoptera results
                    lep_res = res.loc[res.Order == "Lepidoptera"]

                    # Number of ASVs
                    foundasvs = res.shape[0]
                    # Number of OTUs
                    otus = len(res["cluster"].unique())
                    # Number of Lepidoptera OTUs
                    otus_lepidoptera = len(lep_res["cluster"].unique())
                    # Number of BOLD_bins
                    bins = len(res["BOLD_bin"].unique())
                    # Number of Lepidoptera BOLD_bins
                    bins_lep = len(lep_res["BOLD_bin"].unique())
                    # Calculate zero species OTUs
                    otu_n_species = len(
                        res.loc[~res.Species.str.contains("unclassified")]["cluster"].unique()
                    )
                    zero_species_otus = otus - otu_n_species
                    # Calculate zero species OTUs lepidoptera
                    otu_n_species_lep = len(
                        lep_res.loc[~lep_res.Species.str.contains("unclassified")][
                            "cluster"
                        ].unique()
                    )
                    zero_species_otus_lepidoptera = otus_lepidoptera - otu_n_species_lep
                    # Calculate zero bin OTUs
                    otu_n_bin = len(
                        res.loc[~res.BOLD_bin.str.contains("unclassified")]["cluster"].unique()
                    )
                    zero_bin_otus = otus - otu_n_bin
                    # Calculate zero bin OTUs Lepidoptera
                    otu_n_bin_lep = len(
                        lep_res.loc[~lep_res.BOLD_bin.str.contains("unclassified")][
                            "cluster"
                        ].unique()
                    )
                    zero_bin_otus_lepidoptera = otus_lepidoptera - otu_n_bin_lep
                    # Species per OTU
                    species_per_otu = (
                        res.loc[~res.Species.str.contains("unclassified")]
                        .groupby("BOLD_bin")
                        .apply(count_species)
                    )
                    # Lepidoptera species per OTU
                    lepidoptera_species_per_otu = (
                        lep_res.loc[~lep_res.Species.str.contains("unclassified")]
                        .groupby("BOLD_bin")
                        .apply(count_species)
                    )
                    # OTUs per species
                    otus_per_species = (
                        res.loc[~res.Species.str.contains("unclassified")]
                        .groupby("Species")
                        .apply(count_otus)
                    )
                    # Lepidoptera OTUs per species
                    lepidoptera_otus_per_species = (
                        lep_res.loc[~lep_res.Species.str.contains("unclassified")]
                        .groupby("Species")
                        .apply(count_otus)
                    )
                    # Total species
                    species = len(res["Species"].unique())
                    # Total reads
                    reads = res.Size.sum()
                    # Multi OTU species
                    multi_otu_species = otus_per_species.loc[otus_per_species > 1].shape[0]
                    multi_otu_species_lepidoptera = lepidoptera_otus_per_species.loc[
                        lepidoptera_otus_per_species > 1
                    ].shape[0]
                    # FINBOL ASVs
                    fin_asvs = fin_res.shape[0]
                    # FINBOL OTUs
                    fin_otus = len(fin_res["cluster"].unique())
                    # FINBOL BOLD_bins
                    fin_bins = len(fin_res["BOLD_bin"].unique())
                    # FINBOL Species
                    fin_species = len(fin_res["Species"].unique())
                    # FINBOL reads
                    fin_reads = fin_res.Size.sum()
                    out = [
                        h,
                        foundasvs,
                        otus,
                        bins,
                        species,
                        reads,
                        fin_asvs,
                        fin_otus,
                        fin_bins,
                        fin_species,
                        fin_reads,
                        zero_bin_otus,
                        zero_species_otus,
                        multi_otu_species,
                        zero_bin_otus_lepidoptera,
                        zero_species_otus_lepidoptera,
                        multi_otu_species_lepidoptera,
                    ]
                    fhout.write("\t".join([str(x) for x in out]) + "\n")


def main(args):
    sys.stderr.write(f"#Reading 'special' ASVs from {args.special_asv_seqs_file}\n")
    asvid_isspecial = get_special_asvs(args.special_asv_seqs_file)
    sys.stderr.write(f"#Read {len(asvid_isspecial.keys())} ASVs\n")
    sys.stderr.write(f"#Reading uchime results from {args.uchime_file}\n")
    uchime_res = read_uchime(args.uchime_file, nonchims=list(asvid_isspecial.keys()))
    sys.stderr.write(f"#Reading taxonomy from {args.taxonomy_file}\n")
    asv_taxonomy = get_taxonomy(args.taxonomy_file)
    sys.stderr.write("#Evaluating results\n")
    eval_uchime(uchime_res, asv_taxonomy, args.minh)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "special_asv_seqs_file",
        type=str,
        help="Fasta file with 'trusted' ASV seqs (nonchimeric)",
    )
    parser.add_argument(
        "taxonomy_file", type=str, help="File with taxonomic info for ASVs"
    )
    parser.add_argument("uchime_file", type=str, help="Uchime output file")
    parser.add_argument(
        "--minh",
        type=str,
        default="0.001,0.005,0.01,0.02,0.04,0.08,0.16,0.28,0.5,1.0,1000",
        help="Comma separated values of minh to evaluate",
    )
    parser.add_argument(
        "--mindiffs",
        type=str,
        default="0,1,2,3,4",
        help="Comma separated values of mindiffs to evaluate",
    )
    parser.add_argument(
        "--mindiv",
        type=str,
        default="0.5,0.8,1.2,1.5",
        help="Comma separated values of mindiv to evaluate",
    )
    args = parser.parse_args()
    main(args)
