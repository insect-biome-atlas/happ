from snakemake.utils import validate

try:
    validate(config, schema="../schemas/config.schema.yaml", set_default=True)
except Exception as e:
    err = e.args[0]
    items = err.split("\n")
    sys.exit(f"""
    ERROR in config file:
    {items[-2]} {items[1]}

    LIKELY CAUSE:
    This error is likely due to a missing or incorrect value in the config file.
    Make sure that no numeric values are quoted and that lists are set to [] when empty.

    WHAT TO DO:
    For example, to set a list parameter to empty, make sure to set it to [] in 
    the config file instead of just commenting out the list items, e.g.:

    This will raise an error
    taxtools:
        #- "sintax"
        #- "vsearch"
        #- "epa-ng"
        #- "sintax+epa-ng"

    Do this instead:
    taxtools: []
    """)

rule taxonomy_filter:
    output:
        taxonomy="results/neeat/taxonomy_filter/{assignment_rank}/cluster_taxonomy.tsv",
    input:
        taxonomy=config["neeat"]["taxfile"]
    log:
        "logs/neeat/taxonomy_filter/{assignment_rank}.log"
    params:
        src=workflow.source_path("../scripts/neeat/taxonomy_filter.py")
    shell:
        """
        python {params.src} -t {input.taxonomy} -r {wildcards.assignment_rank} -o {output.taxonomy} >{log} 2>&1
        """

def get_taxonomy(wildcards):
    if config["noise_filtering"]["assignment_rank"] == "":
        return config["neeat"]["taxfile"]
    else:
        return expand("results/neeat/taxonomy_filter/{assignment_rank}/cluster_taxonomy.tsv", assignment_rank=config["noise_filtering"]["assignment_rank"])

rule generate_counts_files:
    """
    Generate the counts files for use with neeat filtering
    """
    output:
        cluster_counts="results/neeat/counts/cluster_counts.tsv",
    input:
        counts=config["neeat"]["countsfile"],
        taxonomy=get_taxonomy,
    log:
        "logs/neeat/generate_counts_files.log"
    params:
        meta=config["metadata"]["metadata_file"] if os.path.exists(config["metadata"]["metadata_file"]) else False,
        sample_id_col=config["metadata"]["sample_id_col"],
        sample_type_col=config["metadata"]["sample_type_col"],
        sample_val=config["metadata"]["sample_val"],
    script: 
        "../scripts/neeat/generate_counts_files.R"

checkpoint generate_taxa_seqs:
    """
    Outputs a fasta file for each taxon (e.g. order) in the dataset. Skips taxa with 
    fewer than 2 ASVs.
    """
    output:
        directory("results/neeat/{noise_rank}/fasta"),
        touch("results/neeat/{noise_rank}/single_otus.tsv")
    input:
        taxonomy=get_taxonomy,
        fasta=config["neeat"]["fastafile"],
        counts=rules.generate_counts_files.output.cluster_counts
    log:
        "logs/neeat/generate_taxa_seqs/{noise_rank}.log"
    params:
        src=workflow.source_path("../scripts/neeat/generate_taxa_seqs.py")
    shell:
        """
        python {params.src} -t {input.taxonomy} -f {input.fasta} -c {input.counts} -r {wildcards.noise_rank} -o {output[0]} >{log} 2>&1
        """

#TODO: Make it clear that the cluster_reps.fasta file needs to have headers as:
# '>ASV_id cluster-id'
rule matchlist_vsearch:
    """
    Create a matchlist for sequences in an order using vsearch
    """
    output:
        "results/neeat/{noise_rank}/vsearch/{tax}.matchlist.tsv"
    input:
        "results/neeat/{noise_rank}/fasta/{tax}.fasta"
    log:
        "logs/neeat/matchlist_vsearch/{noise_rank}/{tax}.log"
    conda: config["vsearch-env"]
    container: "docker://quay.io/biocontainers/vsearch:2.29.1--h6a68c12_0"
    params:
        tmpdir = "$TMPDIR/{noise_rank}/{tax}/",
        maxhits = config["noise_filtering"]["max_target_seqs"]
    threads: 4
    wrapper:
        "file:workflow/wrappers/matchlist_vsearch"

rule generate_datasets:
    """
    Output taxonomy and counts files for a certain taxon
    """
    output:
        taxonomy="results/neeat/{noise_rank}/data/{tax}_taxonomy.tsv",
        counts="results/neeat/{noise_rank}/data/{tax}_counts.tsv",
    input:
        counts=rules.generate_counts_files.output.cluster_counts,
        taxonomy=get_taxonomy,
    log:
        "logs/neeat/generate_datasets/{noise_rank}/{tax}.log"
    params:
        split_rank = config["noise_filtering"]["split_rank"]
    threads: 1
    script:
        "../scripts/neeat/generate_datasets.R"

rule generate_aa_seqs:
    """
    Translate nucleotide sequences to amino acid sequences
    """
    output:
        faa="results/neeat/{noise_rank}/faa/{tax}.faa"
    input:
        fasta="results/neeat/{noise_rank}/fasta/{tax}.fasta"
    log:
        "logs/neeat/generate_aa_seqs/{noise_rank}/{tax}.log"
    params:
        codon_table = config["noise_filtering"]["codon_table"],
        codon_start = config["noise_filtering"]["codon_start"]
    script:
        "../scripts/neeat/generate_aa_seqs.R"

rule mafft_align:
    """
    Align protein sequences using MAFFT
    """
    output:
        "results/neeat/{noise_rank}/mafft/{tax}.aligned.faa"
    input:
        rules.generate_aa_seqs.output.faa,
    log:
        "logs/neeat/mafft_align/{noise_rank}/{tax}.log",
    conda:  config["mafft-env"]
    container: "docker://quay.io/biocontainers/mafft:7.525--h031d066_0"
    threads: 4
    wrapper:
        "file:workflow/wrappers/mafft_align"

rule trim_align:
    output:
        nuc="results/neeat/{noise_rank}/trimmed/{tax}.aligned.fasta"
    input:
        nuc="results/neeat/{noise_rank}/fasta/{tax}.fasta"
    log:
        "logs/neeat/trim_align/{noise_rank}/{tax}.log"
    params:
        codon_start = config["noise_filtering"]["codon_start"],
    wrapper:
        "file:workflow/wrappers/trim_align"

rule pal2nal:
    """
    Generate the corresponding nucleotide alignments with pal2nal
    """
    output:
        ensure("results/neeat/{noise_rank}/pal2nal/{tax}.fasta", non_empty=True)
    input:
        pep=rules.mafft_align.output[0],
        nuc=rules.trim_align.output.nuc
    log:
        "logs/neeat/pal2nal/{noise_rank}/{tax}.log"
    conda: config["pal2nal-env"]
    container: "docker://biocontainers/pal2nal:v14.1-2-deb_cv1"
    params:
        codon_table = config["noise_filtering"]["codon_table"],
    wrapper:
        "file:workflow/wrappers/pal2nal"

rule generate_evodistlists:
    """
    generates evolutionary distance files for the evo_filter function of neeat
    """
    output:
        tsv="results/neeat/{noise_rank}/evodist/{tax}_evodistlist.tsv"
    input:
        matchlist=rules.matchlist_vsearch.output[0],
        taxonomy=rules.generate_datasets.output.taxonomy,
        fasta=rules.pal2nal.output[0]
    log:
        "logs/neeat/generate_evodistlists/{noise_rank}/{tax}.log"
    params:
        codon_model = workflow.source_path("../scripts/neeat/codon_model.R"),
    container: "docker://quay.io/biocontainers/r-seqinr:3.4_5--r3.4.1_0"
    conda: config["seqinr-env"]
    threads: 1
    script:
        "../scripts/neeat/generate_evo_dists.R"

def get_checkpoint_files(wildcards):
    checkpoint_dir = checkpoints.generate_taxa_seqs.get(**wildcards).output[0]
    files = expand("results/neeat/{noise_rank}/evodist/{tax}_evodistlist.tsv",
        noise_rank=wildcards.noise_rank, tax=glob_wildcards(os.path.join(checkpoint_dir, "{tax}.fasta")).tax)
    return files

rule collect_files:
    input:
        get_checkpoint_files
    output:
        touch("results/neeat/{noise_rank}/evodist.txt")