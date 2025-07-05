from snakemake.utils import validate
import os

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

for key in [x for x in config.keys() if x.endswith("-env")]:
    pixi_dir = f".pixi/envs/{key.replace('-env', '')}"
    if os.path.exists(pixi_dir) and os.path.isdir(pixi_dir):
        config[key] = pixi_dir

rule all:
    input:
        expand("results/neeat/{noise_rank}/{f}",
            noise_rank = config["noise_filtering"]["split_rank"],
            f = ["noise_filtered_cluster_counts.tsv","noise_filtered_cluster_taxonomy.tsv","discarded_cluster_taxonomy.tsv"]
        )

rule taxonomy_filter:
    message: "Removing clusters unassigned at {wildcards.assignment_rank}"
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
    message: "Generating countsfile for NEEAT"
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
    message: "Generating fasta files for taxa at {wildcards.noise_rank}"
    output:
        directory("results/neeat/{noise_rank}/fasta"),
        touch("results/neeat/{noise_rank}/fasta/singles.tsv")
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

rule matchlist_vsearch:
    """
    Create a matchlist for sequences in an order using vsearch
    """
    message: "Comparing sequences for {wildcards.tax}"
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
    shell:
        """
        vsearch --usearch_global {input} --db {input} --self --id .84 --iddef 1 \
            --userout {output} -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits {params.maxhits} --threads {threads} > {log} 2>&1
        """

rule generate_datasets:
    """
    Output taxonomy and counts files for a certain taxon
    """
    message: "Generating files for {wildcards.tax}"
    output:
        taxonomy="results/neeat/{noise_rank}/data/{tax}_taxonomy.tsv",
        counts="results/neeat/{noise_rank}/data/{tax}_counts.tsv",
    input:
        counts=rules.generate_counts_files.output.cluster_counts,
        taxonomy=get_taxonomy,
    log:
        "logs/neeat/generate_datasets/{noise_rank}/{tax}.log"
    params:
        split_rank = config["noise_filtering"]["split_rank"],
        src=workflow.source_path("../scripts/neeat/generate_datasets.py"),
        outdir=lambda wildcards, output: os.path.dirname(output.taxonomy)
    threads: 1
    shell:
        """
        python {params.src} {input.taxonomy} {input.counts} \
            -r {wildcards.noise_rank} -t {wildcards.tax} -o {params.outdir}
        """

rule generate_aa_seqs:
    """
    Translate nucleotide sequences to amino acid sequences
    """
    message: "Translating sequences for {wildcards.tax} using codon table {params.codon_table}"
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
    message: "Aligning translated sequences with MAFFT for {wildcards.tax}"
    output:
        "results/neeat/{noise_rank}/mafft/{tax}.aligned.faa"
    input:
        rules.generate_aa_seqs.output.faa,
    log:
        "logs/neeat/mafft_align/{noise_rank}/{tax}.log",
    conda:  config["mafft-env"]
    container: "docker://quay.io/biocontainers/mafft:7.525--h031d066_0"
    threads: 4
    shell:
        """
        mafft --auto --thread {threads} {input} > {output} 2>{log}
        """

rule trim_align:
    message: "Trimming alignments for {wildcards.tax}"
    output:
        nuc="results/neeat/{noise_rank}/trimmed/{tax}.aligned.fasta"
    input:
        nuc="results/neeat/{noise_rank}/fasta/{tax}.fasta"
    log:
        "logs/neeat/trim_align/{noise_rank}/{tax}.log"
    params:
        codon_start = config["noise_filtering"]["codon_start"],
    shell:
        """
        seqkit subseq --region {params.codon_start}:-1 {input.nuc} > {output.nuc} 2>{log}
        """

rule pal2nal:
    """
    Generate the corresponding nucleotide alignments with pal2nal
    """
    message: "Generating nucleotide alignments with pal2nal for {wildcards.tax}"
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
    shell:
        """
        pal2nal.pl {input.pep} {input.nuc} -output fasta -codontable {params.codon_table} > {output} 2>{log}
        """

rule generate_evodistlists:
    """
    Generates evolutionary distance files for the evo_filter function of neeat
    """
    message: "Generating evolutionary distances for {wildcards.tax}"
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

def agg_evodist(wc):
    """
    Aggregate the evodistlists for all taxa
    """
    checkpoint_dir = checkpoints.generate_taxa_seqs.get(**wc).output[0]
    return expand("results/neeat/{noise_rank}/evodist/{tax}_evodistlist.tsv",
        noise_rank=wc.noise_rank, tax=glob_wildcards(os.path.join(checkpoint_dir, "{tax}.fasta")).tax)

rule aggregate_evodist:
    input:
        agg_evodist
    output:
        touch("results/neeat/{noise_rank}/evodist.done")

rule generate_neeat_filtered:
    message: "NEEAT filtering sequences for {wildcards.tax}"
    output:
        retained="results/neeat/{noise_rank}/filtered/{tax}_retained.tsv",
        discarded="results/neeat/{noise_rank}/filtered/{tax}_discarded.tsv",
        counts="results/neeat/{noise_rank}/filtered/{tax}_counts.tsv",
    input:
        counts=rules.generate_datasets.output.counts,
        distlist=rules.generate_evodistlists.output[0],
        taxonomy=rules.generate_datasets.output.taxonomy,
    log:
        "logs/neeat/neeat_filter/{noise_rank}/{tax}.log"
    params:
        neeat_filter=workflow.source_path("../scripts/neeat/neeat_filter.R"),
        echo_filter=workflow.source_path("../scripts/neeat/echo_filter.R"),
        evo_filter=workflow.source_path("../scripts/neeat/evo_filter.R"),
        abundance_filter=workflow.source_path("../scripts/neeat/abundance_filter.R"),
        min_match=config["noise_filtering"]["min_match"],
        n_closest=config["noise_filtering"]["n_closest"],
        echo_min_overlap=config["noise_filtering"]["echo_min_overlap"],
        echo_read_ratio_type=config["noise_filtering"]["echo_read_ratio_type"],
        echo_max_read_ratio=config["noise_filtering"]["echo_max_read_ratio"],
        echo_require_corr=config["noise_filtering"]["echo_require_corr"],
        evo_local_min_overlap=config["noise_filtering"]["evo_local_min_overlap"],
        dist_type_local=config["noise_filtering"]["dist_type_local"],
        dist_threshold_local=config["noise_filtering"]["dist_threshold_local"],
        dist_threshold_global=config["noise_filtering"]["dist_threshold_global"],
        abundance_cutoff_type=config["noise_filtering"]["abundance_cutoff_type"],
        abundance_cutoff=config["noise_filtering"]["abundance_cutoff"],
        assignment_rank=config["noise_filtering"]["assignment_rank"]
    script:
        "../scripts/neeat/generate_neeat_filtered.R"

def aggregate_neeat(wc):
    """
    Aggregate the evodistlists for all taxa
    """
    checkpoint_dir = checkpoints.generate_taxa_seqs.get(**wc).output[0]
    retained = expand("results/neeat/{noise_rank}/filtered/{tax}_retained.tsv",
        noise_rank=wc.noise_rank, tax=glob_wildcards(os.path.join(checkpoint_dir, "{tax}.fasta")).tax)
    return retained

rule neeat:
    message: "Aggregating NEEAT filtered files"
    output:
        counts = "results/neeat/{noise_rank}/noise_filtered_cluster_counts.tsv",
        retained = "results/neeat/{noise_rank}/noise_filtered_cluster_taxonomy.tsv",
        discarded = "results/neeat/{noise_rank}/discarded_cluster_taxonomy.tsv"
    input:
        retained=aggregate_neeat,
        counts = config["neeat"]["countsfile"],
        taxonomy=config["neeat"]["taxfile"],
        singles = "results/neeat/{noise_rank}/fasta/singles.tsv"
    log:
        "logs/neeat/{noise_rank}/neeat.log"
    params:
        src=workflow.source_path("../scripts/neeat/neeat.py"),
        outdir = lambda wildcards, output: os.path.dirname(output.counts)
    shell:
        """
        python {params.src} -r {input.retained} -c {input.counts} -t {input.taxonomy} -s {input.singles} -o {params.outdir} > {log} 2>&1
        """
