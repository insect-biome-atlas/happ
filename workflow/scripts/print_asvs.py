#/usr/bin/env python

from argparse import ArgumentParser
from Bio.SeqIO import parse
import sys
from collections import defaultdict

def get_special_asvs(special_asv_seqs_file):
    asvid_isspecial = {}
    for record in parse(special_asv_seqs_file, 'fasta'):
        asvid_isspecial[record.id] = 1
    return asvid_isspecial

def get_taxonomy(taxonomy_file):
    asvid_otu = {}
    asvid_bin = {}
    asvid_species = {}
    asvid_order = {}
    with open(taxonomy_file, 'r') as INFILE:
        for i, line in enumerate(INFILE):
            line = line.rstrip()
            if i == 0:
                continue
            fields = line.rsplit()
            asvid_otu[fields[0]] = fields[1]
            asvid_bin[fields[0]] = fields[10]
            asvid_species[fields[0]] = fields[9]
            asvid_order[fields[0]] = fields[6]
    return asvid_order, asvid_bin, asvid_species, asvid_otu

def get_uchime(uchime_file):
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
    asvid_minh = {}
    asvid_count = {}
    asvid_sumL = {}
    asvid_sumR = {}
    asvid_count = {}
    asvid_div = {}
    with open(uchime_file, 'r') as INFILE:
        for line in INFILE:
            line = line.rstrip()
            fields = line.rsplit()
            asvid = fields[1].split(";")[0]
            size = int(fields[1].split("=")[-1])
            asvid_minh[asvid] = float(fields[0])
            asvid_count[asvid] = size
            asvid_sumL[asvid] = int(fields[10]) + int(fields[11]) + int(fields[12])
            asvid_sumR[asvid] = int(fields[13]) + int(fields[14]) + int(fields[15])
            asvid_div[asvid] = float(fields[16])
    return asvid_div, asvid_sumR, asvid_sumL, asvid_count, asvid_minh


def eval_uchime(asvid_isspecial, asvid_order, asvid_bin, asvid_species, asvid_otu,
                asvid_div, asvid_sumR, asvid_sumL, asvid_count, asvid_minh,
                levels=[0.001, 0.005, 0.01, 0.02, 0.04, 0.08, 0.16, 0.28, 0.5, 1.0, 1000],):
    with sys.stdout as fhout:
        fhout.write("minh\tASVs\totus\tbins\tspecies\treads\tFIN_ASVs\tFIN_otus\tFIN_bins\tFIN_species\tFIN_reads\tzero_bin_otus\tzero_species_otus\tmulti_otu_species\tzero_bin_otus_lepidoptera\tzero_species_otus_lepidoptera\tmulti_otu_species_lepidoptera\n")
        for level in levels:
            foundasvs = {}; foundotus = {}; foundbins = {}; foundspecies = {}; foundreads = 0
            foundasvs_spec = {}; foundotus_spec = {}; foundbins_spec = {}; foundspecies_spec = {}; foundreads_spec = 0

            foundotus_lepidoptera = {}
            otu_n_species = defaultdict(lambda: defaultdict(dict))
            species_n_otu = defaultdict(lambda: defaultdict(dict))
            otu_n_species_lepidoptera = defaultdict(lambda: defaultdict(dict))
            species_n_otu_lepidoptera = defaultdict(lambda: defaultdict(dict))
            otu_n_bin = defaultdict(lambda: defaultdict(dict))
            bin_n_otu = defaultdict(lambda: defaultdict(dict))
            otu_n_bin_lepidoptera = defaultdict(lambda: defaultdict(dict))
            bin_n_otu_lepidoptera = {}

            for asvid in asvid_otu:
                try:
                    asvid_minh[asvid]
                except KeyError:
                    continue
                if asvid_minh[asvid] < level:
                    foundasvs[asvid] = 1
                    foundotus[asvid_otu[asvid]] = 1
                    foundbins[asvid_bin[asvid]] = 1
                    foundspecies[asvid_species[asvid]] = 1
                    foundreads+=asvid_count[asvid]
                    if "unclassified" not in asvid_species[asvid]:
                        otu_n_species[asvid_otu[asvid]][asvid_species[asvid]] = 1
                        species_n_otu[asvid_species[asvid]][asvid_otu[asvid]] = 1

                        if asvid_order[asvid] == "Lepidoptera":
                            otu_n_species_lepidoptera[asvid_otu[asvid]][asvid_species[asvid]] = 1
                            species_n_otu_lepidoptera[asvid_species[asvid]][asvid_otu[asvid]] = 1
                    if "unclassified" not in asvid_bin[asvid]:
                        otu_n_bin[asvid_otu[asvid]][asvid_bin[asvid]] = 1
                        bin_n_otu[asvid_bin[asvid]][asvid_otu[asvid]] = 1
                        if asvid_order[asvid] == "Lepidoptera":
                            otu_n_bin_lepidoptera[asvid_otu[asvid]][asvid_bin[asvid]] = 1
                            bin_n_otu_lepidoptera[asvid_bin[asvid]][asvid_otu[asvid]] = 1
                    if asvid_order[asvid] == "Lepidoptera":
                        foundotus_lepidoptera[asvid_otu[asvid]] = 1
                    if asvid in asvid_isspecial.keys():
                        foundasvs_spec[asvid] = 1
                        foundotus_spec[asvid_otu[asvid]] = 1
                        foundbins_spec[asvid_bin[asvid]] = 1
                        foundspecies_spec[asvid_species[asvid]] = 1
                        foundreads_spec+=asvid_count[asvid]
            stats = []
            stats.append(len(foundasvs))
            stats.append(len(foundotus))
            stats.append(len(foundbins))
            stats.append(len(foundspecies))
            stats.append(foundreads)
            stats.append(len(foundasvs_spec))
            stats.append(len(foundotus_spec))
            stats.append(len(foundbins_spec))
            stats.append(len(foundspecies_spec))
            stats.append(foundreads_spec)
            stats.append(len(foundotus)-len(otu_n_bin)) # zero_bin_otu
            stats.append(len(foundotus)-len(otu_n_species)) # zero_species_otu
            multi_otu_species = 0
            for key in species_n_otu.keys():
                if len(species_n_otu[key]) > 1:
                    multi_otu_species+=1
            stats.append(multi_otu_species)
            stats.append(len(foundotus_lepidoptera) - len(otu_n_bin_lepidoptera)) # zero_bin_otu lepidoptera
            stats.append(len(foundotus_lepidoptera) - len(otu_n_species_lepidoptera)) # zero_species_otu lepidoptera
            multi_otu_species_lepidoptera = 0
            for key in species_n_otu_lepidoptera.keys():
                if len(species_n_otu_lepidoptera[key]) > 1:
                    multi_otu_species_lepidoptera+=1
            stats.append(multi_otu_species_lepidoptera)
            fhout.write([level]+"\t".join(stats)+"\n")


def main(args):
    asvid_isspecial = get_special_asvs(args.special_asv_seqs_file)
    asvid_order, asvid_bin, asvid_species, asvid_otu = get_taxonomy(args.taxonomy_file)
    asvid_div, asvid_sumR, asvid_sumL, asvid_count, asvid_minh = get_uchime(args.uchime_fike)
    eval_uchime(asvid_isspecial, asvid_order, asvid_bin, asvid_species, asvid_otu,
                asvid_div, asvid_sumR, asvid_sumL, asvid_count, asvid_minh)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("special_asv_seqs_file", type=str,
                      help="Fasta file with 'trusted' ASV seqs (nonchimeric)")
    parser.add_argument("taxonomy_file", type=str,
                      help="File with taxonomic info for ASVs")
    parser.add_argument("uchime_file", type=str,
                      help="Uchime output file")
    args = parser.parse_args()
    main(args)
