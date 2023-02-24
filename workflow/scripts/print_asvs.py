### files ###

special_asv_seqs_file = "iba_data/AllSamples_2209.221216.opticlust.025/asv_seqs_in_finbol.fasta"
taxonomy_file = "iba_data/AllSamples_2209.221216.opticlust.025/asv_reps.taxonomy.tsv"
uchime_file = "iba_data/AllSamples_2209.221216.opticlust.025/uchimeout.txt"
#samplewise_chimera_file = "iba_data/AllSamples_2209.221216.opticlust.025/samplewise_chimeras.tsv"
#all_asv_seqs_file = "iba_data/AllSamples_2209.221216.opticlust.025/asv_seqs.fasta"
#aligned_asvs_file = "iba_data/AllSamples_2209.221216.opticlust.025/trouble_shooting/Plutella_xylostella_ASVs_aln_mafft.fasta"

### parameters ###

#query_taxonomy = "Plutella xylostella"
#query_taxonomy = "Phaonia lugubris"
#query_taxonomy_level = 9
#minh_cutoff = 0.28
#minh_cutoff = 0.005

### main code ###

def get_special_asvs():
    pass

def get_taxonomy():
    pass

def get_uchime():
    pass

def eval_uchime():
    levels = [0.001, 0.005, 0.01, 0.02, 0.04, 0.08, 0.16, 0.28, 0.50, 1.0, 1000] # Minimum score (h) to be classified as chimera.
    print("minh\tASVs\totus\tbins\tspecies\treads\tFIN_ASVs\tFIN_otus\tFIN_bins\tFIN_species\tFIN_reads\tzero_bin_otus\tzero_species_otus\tmulti_otu_species\tzero_bin_otus_lepidoptera\tzero_species_otus_lepidoptera\tmulti_otu_species_lepidoptera")
    for level in levels:
        foundasvs = {}
        foundotus = {}
        foundbins = {}
        foundspecies = {}
        foundreads = 0

        foundasvs_spec = {}
        foundotus_spec = {}
        foundbins_spec = {}
        foundspecies_spec = {}
        foundreads_spec = 0

        foundotus_lepidoptera = {}
        otu_n_species = {}
        species_n_otu = {}
        otu_n_species_lepidoptera = {}
        species_n_otu_lepidoptera = {}
        otu_n_bin = {}
        bin_n_otu = {}
        otu_n_bin_lepidoptera = {}
        bin_n_otu_lepidoptera = {}

        for asvid in asvid_otu:
            #if asvid_order[asvid] != "Lepidoptera":
            #    continue
            if asvid in asvid_minh:
                if asvid_minh[asvid] < level:
                    foundasvs[asvid] = 1
                    foundotus[asvid_otu[asvid]] = 1
                    foundbins[asvid_bin[asvid]] = 1
                    foundspecies[asvid_species[asvid]] = 1
                    foundreads += asvid_count[asvid]
                    if "_X" not in asvid_species[asvid] and "unclassified" not in asvid_species[asvid]:
                        otu_n_species.setdefault(asvid_otu[asvid], {}
