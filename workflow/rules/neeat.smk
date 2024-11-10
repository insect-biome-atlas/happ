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
    run:
        import pandas as pd
        assignment_rank = wildcards.assignment_rank
        taxdf = pd.read_csv(input.taxonomy, sep="\t", index_col=0)
        taxdf = taxdf.loc[(~taxdf[assignment_rank].str.contains("_X+$"))&(~taxdf[assignment_rank].str.startswith("unclassified"))]
        taxdf.to_csv(output.taxonomy, sep="\t")

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
    run:
        import subprocess
        import pandas as pd
        rank = wildcards.noise_rank
        outdir=output[0]
        os.makedirs(outdir, exist_ok=True)
        # read in the taxonomy file
        taxdf = pd.read_csv(input.taxonomy[0], sep="\t", index_col=0)
        # extract representative ASVs
        taxdf = taxdf.loc[taxdf["representative"]==1]
        # also make sure that only clusters that remain in the counts file are used
        clusters_w_counts = []
        with open(input.counts, 'r') as fhin:
            for line in fhin:
                clusters_w_counts.append(line.split("\t")[0])
        taxdf = taxdf.loc[taxdf["cluster"].isin(clusters_w_counts)]
        # taxa is the unique set of values for split_rank
        taxa = taxdf[rank].unique()
        singles = pd.DataFrame()
        # iterate over taxa
        for tax in taxa:
            # get asvs and cluster designation
            asvs = taxdf.loc[taxdf[rank]==tax, "cluster"]
            # skip taxa with fewer than 2 ASVs
            if len(asvs) < 2:
                singles = pd.concat([singles, asvs])
                continue
            # generate strings with '<asv> <cluster>' format to match fasta headers
            asvs = [f"{x[0]} {x[1]}" for x in list(zip(asvs.index, asvs))]
            # write asvs to file
            with open(f"{output[0]}/{tax}.ids", "w") as f:
                _ = f.write("\n".join(asvs))
            # use seqkit to extract sequences
            with open(f"{outdir}/{tax}.fasta", 'w') as fhout:
                cmd = ["seqkit", "grep", "-f",f"{output[0]}/{tax}.ids","-n", "--quiet", input.fasta]
                _ = subprocess.run(cmd, stdout=fhout)
                os.remove(f"{output[0]}/{tax}.ids")
        singles.to_csv(output[1], sep="\t")



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
    resources:
        tasks = 4
    shell:
        """
        mkdir -p {params.tmpdir}
        cat {input} | sed 's/>.\+ />/g' > {params.tmpdir}/cluster_reps.fasta
        vsearch --usearch_global {params.tmpdir}/cluster_reps.fasta --db {params.tmpdir}/cluster_reps.fasta --self --id .84 --iddef 1 \
            --userout {output} -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits {params.maxhits} --threads {resources.tasks} > {log} 2>&1
        rm -rf {params.tmpdir}
        """

rule generate_datasets:
    """
    Output taxonomy and counts files for a certain taxon
    """
    output:
        taxonomy="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/data/{tax}_taxonomy.tsv",
        counts="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/data/{tax}_counts.tsv",
        #cal_counts="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/data/{tax}_cal_counts.tsv",
        #tot_prop_counts="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/data/{tax}_tot_prop_counts.tsv",
        #sample_prop_counts="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/data/{tax}_sample_prop_counts.tsv",
    input:
        counts=rules.generate_counts_files.output.cluster_counts,
        #cal_counts=rules.generate_counts_files.output.calibrated_counts,
        #tot_prop_counts=rules.generate_counts_files.output.tot_proportional_counts,
        #sample_prop_counts=rules.generate_counts_files.output.sample_proportional_counts,
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
    resources:
        tasks = 4,
        cpus_per_task = 1
    shell:
        "mafft --auto --thread {resources.tasks} {input} > {output} 2>{log}"

rule trim_align:
    output:
        nuc="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/trimmed/{tax}.aligned.faa"
    input:
        nuc="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/fasta/{tax}.fasta"
    params:
        codon_start = config["noise_filtering"]["codon_start"],
    shell:
        """
        seqkit subseq --region {params.codon_start}:-1 {input.nuc} > {output.nuc}
        """

rule pal2nal:
    """
    Generate the corresponding nucleotide alignments with pal2nal
    """
    output:
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/pal2nal/{tax}.fasta"
    input:
        pep=rules.mafft_align.output[0],
        nuc=rules.trim_align.output.nuc
    log:
        "logs/noise_filtering/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{noise_rank}/{tax}.log"
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

def concat(files):
    """
    Concatenates multiple tab-separated value (TSV) files into a single DataFrame.

    Args:
        files (list of str): List of file paths to the TSV files to be concatenated.

    Returns:
        pandas.DataFrame: A DataFrame containing the concatenated data from all input files.

    The function reads each file into a DataFrame, ensuring that all DataFrames have the same columns
    by using the columns from the first file. It then concatenates all DataFrames along the row axis.
    """
    df = pd.DataFrame()
    for i, f in enumerate(files):
            _df = pd.read_csv(f, sep="\t", index_col=0)
            if i==0:
                cols = _df.columns
            else:
                _df = _df.loc[:, cols]
            df = pd.concat([df, _df], axis=0)
    return df

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
    run:
        import pandas as pd
        # filter cluster taxonomy
        taxdf = pd.read_csv(input.taxonomy, sep="\t", index_col=0)
        counts = pd.read_csv(input.counts, sep="\t", index_col=0)
        retained = concat(input.retained)
        singles = pd.read_csv(input.singles, sep="\t", index_col=0, header=0, names=["ASV","cluster"])
        if len(singles) > 0:
            singles = taxdf.loc[singles.index].reset_index().set_index("cluster")
            singles = singles.loc[:, retained.columns]
            retained = pd.concat([retained, singles])
        # write discarded
        taxdf.loc[~taxdf["cluster"].isin(retained.index)].to_csv(output.discarded, sep="\t")
        # write retained
        taxdf.loc[taxdf["cluster"].isin(retained.index)].to_csv(output.retained, sep="\t")
        # filter consensus taxonomy
        cons_tax = pd.read_csv(input.consensus, sep="\t", index_col=0)
        cons_tax.loc[retained.index].to_csv(output.cons_retained, sep="\t")
        # filter counts
        counts.loc[retained.index].to_csv(output.counts, sep="\t")

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