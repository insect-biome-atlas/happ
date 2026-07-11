localrules:
    generate_counts_files,
    generate_taxa_seqs,
    taxonomy_filter,
    trim_align,
    generate_aa_seqs,
    neeat,


rule taxonomy_filter:
    input:
        taxonomy="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_taxonomy.tsv",
    output:
        taxonomy="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/taxonomy_filter/{assignment_rank}/cluster_taxonomy.tsv",
    log:
        "logs/taxonomy_filter/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{assignment_rank}.log",
    params:
        src=workflow.source_path("../scripts/neeat/taxonomy_filter.py"),
    message:
        "Removing clusters unassigned at {wildcards.assignment_rank}"
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
        return expand(
            "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_taxonomy.tsv",
            tool=wildcards.tool,
            rundir=wildcards.rundir,
            chimera_run=wildcards.chimera_run,
            chimdir=wildcards.chimdir,
            rank=wildcards.rank,
            run_name=wildcards.run_name,
        )
    else:
        return expand(
            "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/taxonomy_filter/{assignment_rank}/cluster_taxonomy.tsv",
            tool=wildcards.tool,
            rundir=wildcards.rundir,
            chimera_run=wildcards.chimera_run,
            chimdir=wildcards.chimdir,
            rank=wildcards.rank,
            run_name=wildcards.run_name,
            assignment_rank=config["noise_filtering"]["assignment_rank"],
        )


# Generate the counts files for use with neeat filtering
rule generate_counts_files:
    input:
        counts="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_counts.tsv",
        taxonomy=get_taxonomy,
    output:
        cluster_counts="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/counts/cluster_counts.tsv",
    log:
        "logs/generate_counts_files/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}.log",
    params:
        meta=(
            config["metadata"]["metadata_file"]
            if os.path.exists(config["metadata"]["metadata_file"])
            else False
        ),
        sample_id_col=config["metadata"]["sample_id_col"],
        sample_type_col=config["metadata"]["sample_type_col"],
        sample_val=config["metadata"]["sample_val"],
    message:
        "Generating countsfile for NEEAT"
    script:
        "../scripts/neeat/generate_counts_files.R"


# Outputs a fasta file for each taxon (e.g. order) in the dataset. Skips taxa
# with fewer than 2 ASVs.
checkpoint generate_taxa_seqs:
    input:
        taxonomy=get_taxonomy,
        fasta="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_reps.fasta",
        counts=rules.generate_counts_files.output.cluster_counts,
    output:
        directory(
            "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/fasta"
        ),
        touch(
            "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/fasta/singles.tsv"
        ),
    log:
        "logs/generate_taxa_seqs/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{noise_rank}.log",
    params:
        src=workflow.source_path("../scripts/neeat/generate_taxa_seqs.py"),
    message:
        "Generating fasta files for taxa at {wildcards.noise_rank}"
    shell:
        """
        python {params.src} -t {input.taxonomy} -f {input.fasta} -c {input.counts} -r {wildcards.noise_rank} -o {output[0]} >{log} 2>&1
        """


# Create a matchlist for sequences in an order using vsearch
rule matchlist_vsearch:
    input:
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/fasta/{tax}.fasta",
    output:
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/vsearch/{tax}.matchlist.tsv",
    log:
        "logs/matchlist_vsearch/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{noise_rank}/{tax}.log",
    conda:
        config["vsearch-env"]
    container:
        "docker://quay.io/biocontainers/vsearch:2.29.1--h6a68c12_0"
    threads: 4
    params:
        tmpdir="$TMPDIR/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/{noise_rank}/{tax}/",
        maxhits=config["noise_filtering"]["max_target_seqs"],
    message:
        "Comparing sequences for {wildcards.tax}"
    shell:
        """
        vsearch --usearch_global {input} --db {input} --self --id .84 --iddef 1 \
            --userout {output} -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits {params.maxhits} --threads {threads} > {log} 2>&1
        """


# Output taxonomy and counts files for a certain taxon
rule generate_datasets:
    input:
        counts=rules.generate_counts_files.output.cluster_counts,
        taxonomy=get_taxonomy,
    output:
        taxonomy="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/data/{tax}_taxonomy.tsv",
        counts="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/data/{tax}_counts.tsv",
    log:
        "logs/generate_datasets/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{noise_rank}/{tax}.log",
    threads: 1
    params:
        split_rank=config["noise_filtering"]["split_rank"],
        src=workflow.source_path("../scripts/neeat/generate_datasets.py"),
        outdir=lambda wildcards, output: os.path.dirname(output.taxonomy),
    message:
        "Generating files for {wildcards.tax}"
    shell:
        """
        python {params.src} {input.taxonomy} {input.counts} \
            -r {wildcards.noise_rank} -t {wildcards.tax} -o {params.outdir}
        """


# Translate nucleotide sequences to amino acid sequences
rule generate_aa_seqs:
    input:
        fasta="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/fasta/{tax}.fasta",
    output:
        faa="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/faa/{tax}.faa",
    log:
        "logs/generate_aa_seqs/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{noise_rank}/{tax}.log",
    params:
        codon_table=config["noise_filtering"]["codon_table"],
        codon_start=config["noise_filtering"]["codon_start"],
    message:
        "Translating sequences for {wildcards.tax} using codon table {params.codon_table}"
    script:
        "../scripts/neeat/generate_aa_seqs.R"


# Align protein sequences using MAFFT
rule mafft_align:
    input:
        rules.generate_aa_seqs.output.faa,
    output:
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/mafft/{tax}.aligned.faa",
    log:
        "logs/mafft_align/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{noise_rank}/{tax}.log",
    conda:
        config["mafft-env"]
    container:
        "docker://quay.io/biocontainers/mafft:7.525--h031d066_0"
    threads: 4
    message:
        "Aligning translated sequences with MAFFT for {wildcards.tax}"
    shell:
        """
        mafft --auto --thread {threads} {input} > {output} 2>{log}
        """


rule trim_align:
    input:
        nuc="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/fasta/{tax}.fasta",
    output:
        nuc="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/trimmed/{tax}.aligned.fasta",
    log:
        "logs/trim_align/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{noise_rank}/{tax}.log",
    params:
        codon_start=config["noise_filtering"]["codon_start"],
    message:
        "Trimming alignments for {wildcards.tax}"
    shell:
        """
        seqkit subseq --region {params.codon_start}:-1 {input.nuc} > {output.nuc} 2>{log}
        """


# Generate the corresponding nucleotide alignments with pal2nal
rule pal2nal:
    input:
        pep=rules.mafft_align.output[0],
        nuc=rules.trim_align.output.nuc,
    output:
        ensure(
            "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/pal2nal/{tax}.fasta",
            non_empty=True,
        ),
    log:
        "logs/pal2nal/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{noise_rank}/{tax}.log",
    conda:
        config["pal2nal-env"]
    container:
        "docker://biocontainers/pal2nal:v14.1-2-deb_cv1"
    params:
        codon_table=config["noise_filtering"]["codon_table"],
    message:
        "Generating nucleotide alignments with pal2nal for {wildcards.tax}"
    shell:
        """
        pal2nal.pl {input.pep} {input.nuc} -output fasta -codontable {params.codon_table} > {output} 2>{log}
        """


# generates evolutionary distance files for the evo_filter function of neeat
rule generate_evodistlists:
    input:
        matchlist=rules.matchlist_vsearch.output[0],
        taxonomy=rules.generate_datasets.output.taxonomy,
        fasta=rules.pal2nal.output[0],
    output:
        tsv="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/evodist/{tax}_evodistlist.tsv",
    log:
        "logs/generate_evodistlists/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{noise_rank}/{tax}.log",
    conda:
        config["seqinr-env"]
    container:
        "docker://quay.io/biocontainers/r-seqinr:3.4_5--r3.4.1_0"
    threads: 1
    params:
        codon_model=workflow.source_path("../scripts/neeat/codon_model.R"),
    message:
        "Generating evolutionary distances for {wildcards.tax}"
    script:
        "../scripts/neeat/generate_evo_dists.R"


def agg_evodist(wc):
    """
    Aggregate the evodistlists for all taxa
    """
    checkpoint_dir = checkpoints.generate_taxa_seqs.get(**wc).output[0]
    return expand(
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/evodist/{tax}_evodistlist.tsv",
        tool=wc.tool,
        rundir=wc.rundir,
        chimera_run=wc.chimera_run,
        chimdir=wc.chimdir,
        rank=wc.rank,
        run_name=wc.run_name,
        noise_rank=wc.noise_rank,
        tax=glob_wildcards(os.path.join(checkpoint_dir, "{tax}.fasta")).tax,
    )


rule aggregate_evodist:
    input:
        agg_evodist,
    output:
        touch(
            "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/evodist.done"
        ),


rule generate_neeat_filtered:
    input:
        counts=rules.generate_datasets.output.counts,
        distlist=rules.generate_evodistlists.output[0],
        taxonomy=rules.generate_datasets.output.taxonomy,
    output:
        retained="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/filtered/{tax}_retained.tsv",
        discarded="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/filtered/{tax}_discarded.tsv",
        counts="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/filtered/{tax}_counts.tsv",
    log:
        "logs/neeat_filter/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{noise_rank}/{tax}.log",
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
        assignment_rank=config["noise_filtering"]["assignment_rank"],
    message:
        "NEEAT filtering sequences for {wildcards.tax}"
    script:
        "../scripts/neeat/generate_neeat_filtered.R"


def aggregate_neeat(wc):
    """
    Aggregate the evodistlists for all taxa
    """
    checkpoint_dir = checkpoints.generate_taxa_seqs.get(**wc).output[0]
    retained = expand(
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/filtered/{tax}_retained.tsv",
        tool=wc.tool,
        rundir=wc.rundir,
        chimera_run=wc.chimera_run,
        chimdir=wc.chimdir,
        rank=wc.rank,
        run_name=wc.run_name,
        noise_rank=wc.noise_rank,
        tax=glob_wildcards(os.path.join(checkpoint_dir, "{tax}.fasta")).tax,
    )
    return retained


rule neeat:
    input:
        retained=aggregate_neeat,
        counts="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_counts.tsv",
        taxonomy="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_taxonomy.tsv",
        consensus="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_consensus_taxonomy.tsv",
        singles="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/fasta/singles.tsv",
    output:
        counts="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/noise_filtered_cluster_counts.tsv",
        retained="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/noise_filtered_cluster_taxonomy.tsv",
        cons_retained="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/noise_filtered_cluster_consensus_taxonomy.tsv",
        discarded="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/discarded_cluster_taxonomy.tsv",
    log:
        "logs/neeat/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{noise_rank}.log",
    params:
        src=workflow.source_path("../scripts/neeat/neeat.py"),
        outdir=lambda wildcards, output: os.path.dirname(output.counts),
    message:
        "Aggregating NEEAT filtered files"
    shell:
        """
        python {params.src} -r {input.retained} -c {input.counts} -t {input.taxonomy} --consensus_taxonomy {input.consensus} -s {input.singles} -o {params.outdir} > {log} 2>&1
        """


# Calculate precision and recall for the clusters
rule noise_filtered_precision_recall:
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
        ignore_taxa=config["ignore_taxa"],
        ignore_rank=config["ignore_rank"],
    shell:
        """
        python {params.src} {input[0]} {input[0]} --rank {params.eval_rank} --order_level {output.txt_order} --ignore_taxa {params.ignore_taxa} --ignore_rank {params.ignore_rank} > {output.txt} 2>{log}
        """
