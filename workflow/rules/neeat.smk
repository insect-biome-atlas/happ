localrules:
    generate_counts_files,
    generate_taxa_seqs,
    taxonomy_filter,
    trim_align,
    generate_aa_seqs,
    neeat,


rule taxonomy_filter:
    output:
        taxonomy="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/taxonomy_filter/{assignment_rank}/cluster_taxonomy.tsv",
    input:
        taxonomy="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_taxonomy.tsv"
    log:
        "logs/taxonomy_filter/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{assignment_rank}.log"
    params:
        src=workflow.source_path("../scripts/neeat/taxonomy_filter.py")
    shell:
        """
        python {params.src} -t {input.taxonomy} -r {wildcards.assignment_rank} -o {output.taxonomy} >{log} 2>&1
        """

# This function `get_taxonomy` generates the file paths for taxonomy results based on the provided wildcards and configuration.
# It checks the `assignment_rank` in the configuration under `noise_filtering`.
# If `assignment_rank` is an empty string, it returns the path for the taxonomy file without the `neeat/taxonomy_filter` directory.
# Otherwise, it includes the `neeat/taxonomy_filter` directory and the `assignment_rank` in the path.
# The paths are generated using the `expand` function with the wildcards: tool, rundir, chimera_run, chimdir, rank, and run_name.
def get_taxonomy(wildcards):
    if config["noise_filtering"]["assignment_rank"] == "":
        return expand("results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_taxonomy.tsv",
            tool=wildcards.tool, rundir=wildcards.rundir, chimera_run=wildcards.chimera_run, chimdir=wildcards.chimdir, rank=wildcards.rank, run_name=wildcards.run_name)
    else:
        return expand("results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/taxonomy_filter/{assignment_rank}/cluster_taxonomy.tsv",
            tool=wildcards.tool, rundir=wildcards.rundir, chimera_run=wildcards.chimera_run, chimdir=wildcards.chimdir, rank=wildcards.rank, run_name=wildcards.run_name, assignment_rank=config["noise_filtering"]["assignment_rank"])

rule generate_counts_files:
    """
    Generate the counts files for use with neeat filtering
    """
    output:
        cluster_counts="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/counts/cluster_counts.tsv",
    input:
        counts="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_counts.tsv",
        taxonomy=get_taxonomy,
    log:
        "logs/generate_counts_files/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}.log"
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
        directory("results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/fasta"),
        touch("results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/single_otus.tsv")
    input:
        taxonomy=get_taxonomy,
        fasta="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_reps.fasta",
        counts=rules.generate_counts_files.output.cluster_counts
    log:
        "logs/generate_taxa_seqs/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{noise_rank}.log"
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
    output:
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/vsearch/{tax}.matchlist.tsv"
    input:
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/fasta/{tax}.fasta"
    log:
        "logs/matchlist_vsearch/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{noise_rank}/{tax}.log"
    conda: config["vsearch-env"]
    container: "docker://quay.io/biocontainers/vsearch:2.29.1--h6a68c12_0"
    params:
        tmpdir = "$TMPDIR/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/{noise_rank}/{tax}/",
        maxhits = config["noise_filtering"]["max_target_seqs"]
    threads: 4
    wrapper:
        "file:workflow/wrappers/matchlist_vsearch"

rule generate_datasets:
    """
    Output taxonomy and counts files for a certain taxon
    """
    output:
        taxonomy="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/data/{tax}_taxonomy.tsv",
        counts="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/data/{tax}_counts.tsv",
    input:
        counts=rules.generate_counts_files.output.cluster_counts,
        taxonomy=get_taxonomy,
    log:
        "logs/generate_datasets/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{noise_rank}/{tax}.log"
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
        faa="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/faa/{tax}.faa"
    input:
        fasta="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/fasta/{tax}.fasta"
    log:
        "logs/generate_aa_seqs/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{noise_rank}/{tax}.log"
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
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/mafft/{tax}.aligned.faa"
    input:
        rules.generate_aa_seqs.output.faa,
    log:
        "logs/mafft_align/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{noise_rank}/{tax}.log",
    conda:  config["mafft-env"]
    container: "docker://quay.io/biocontainers/mafft:7.525--h031d066_0"
    threads: 4
    wrapper:
        "file:workflow/wrappers/mafft_align"

rule trim_align:
    output:
        nuc="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/trimmed/{tax}.aligned.fasta"
    input:
        nuc="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/fasta/{tax}.fasta"
    log:
        "logs/trim_align/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{noise_rank}/{tax}.log"
    params:
        codon_start = config["noise_filtering"]["codon_start"],
    wrapper:
        "file:workflow/wrappers/trim_align"

rule pal2nal:
    """
    Generate the corresponding nucleotide alignments with pal2nal
    """
    output:
        ensure("results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/pal2nal/{tax}.fasta", non_empty=True)
    input:
        pep=rules.mafft_align.output[0],
        nuc=rules.trim_align.output.nuc
    log:
        "logs/pal2nal/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{noise_rank}/{tax}.log"
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
        tsv="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/evodist/{tax}_evodistlist.tsv"
    input:
        matchlist=rules.matchlist_vsearch.output[0],
        taxonomy=rules.generate_datasets.output.taxonomy,
        fasta=rules.pal2nal.output[0]
    log:
        "logs/generate_evodistlists/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{noise_rank}/{tax}.log"
    params:
        codon_model = workflow.source_path("../scripts/neeat/codon_model.R"),
    container: "docker://quay.io/biocontainers/r-seqinr:3.4_5--r3.4.1_0"
    conda: config["seqinr-env"]
    threads: 1
    script:
        "../scripts/neeat/generate_evo_dists.R"

def agg_evodist(wc):
    """
    Aggregate the evodistlists for all taxa
    """
    checkpoint_dir = checkpoints.generate_taxa_seqs.get(**wc).output[0]
    return expand("results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/evodist/{tax}_evodistlist.tsv",
        tool=wc.tool, rundir=wc.rundir, chimera_run=wc.chimera_run, chimdir=wc.chimdir, rank=wc.rank,
        run_name=wc.run_name, noise_rank=wc.noise_rank, tax=glob_wildcards(os.path.join(checkpoint_dir, "{tax}.fasta")).tax)

rule aggregate_evodist:
    input:
        agg_evodist
    output:
        touch("results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/evodist.done")

rule generate_neeat_filtered:
    output:
        retained="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/filtered/{tax}_retained.tsv",
        discarded="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/filtered/{tax}_discarded.tsv",
        counts="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/filtered/{tax}_counts.tsv",
    input:
        counts=rules.generate_datasets.output.counts,
        distlist=rules.generate_evodistlists.output[0],
        taxonomy=rules.generate_datasets.output.taxonomy,
    log:
        "logs/neeat_filter/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{noise_rank}/{tax}.log"
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
    retained = expand("results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/filtered/{tax}_retained.tsv",
        tool=wc.tool, rundir=wc.rundir, chimera_run=wc.chimera_run, chimdir=wc.chimdir, rank=wc.rank,
        run_name=wc.run_name, noise_rank=wc.noise_rank, tax=glob_wildcards(os.path.join(checkpoint_dir, "{tax}.fasta")).tax)
    return retained

rule neeat:
    output:
        counts = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/noise_filtered_cluster_counts.tsv",
        retained = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/noise_filtered_cluster_taxonomy.tsv",
        cons_retained = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/noise_filtered_cluster_consensus_taxonomy.tsv",
        discarded = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/discarded_cluster_taxonomy.tsv"
    input:
        retained=aggregate_neeat,
        counts = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_counts.tsv",
        taxonomy="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_taxonomy.tsv",
        consensus="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_consensus_taxonomy.tsv",
        singles = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/single_otus.tsv"
    log:
        "logs/neeat/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{noise_rank}.log"
    params:
        src=workflow.source_path("../scripts/neeat/neeat.py"),
        outdir = lambda wildcards, output: os.path.dirname(output.counts)
    shell:
        """
        python {params.src} -r {input.retained} -c {input.counts} -t {input.taxonomy} --consensus_taxonomy {input.consensus} -s {input.singles} -o {params.outdir} > {log} 2>&1
        """

rule noise_filtered_precision_recall:
    """
    Calculate precision and recall for the clusters
    """
    input:
        rules.neeat.output.retained,
    output:
        txt="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/noise_filtered_precision_recall.txt",
        txt_order="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/noise_filtered_precision_recall.order.txt",
    log:
        "logs/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/precision_recall.log",
    params:
        src=workflow.source_path("../scripts/evaluate_clusters.py"),
        eval_rank=config["evaluation_rank"],
    shell:
        """
        python {params.src} {input[0]} {input[0]} --rank {params.eval_rank} --order_level {output.txt_order} > {output.txt} 2>{log}
        """