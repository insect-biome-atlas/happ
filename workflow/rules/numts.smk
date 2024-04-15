import pandas as pd

for dataset in config["input"].keys():
    if not "spikein_clusters" in config["input"][dataset].keys():
        config["input"][dataset]["spikein_clusters"] = []
    if not "non_numt_ASVs" in config["input"][dataset].keys():
        config["input"][dataset]["non_numt_ASVs"] = ""
    if not "consensus_taxonomy" in config["input"][dataset].keys():
        config["input"][dataset]["consensus_taxonomy"] = ""

def mem_allowed(wildcards, threads):
    return max(threads * 6400, 6400)

localrules: 
    generate_cluster_analysis,
    generate_order_seqs,
    generate_trimmed_seq,
    generate_aa_seqs,
    abundance_filter,
    evaluate_order,
    filtered

def numts_input(wildcards):
    input = []
    for dataset in config["input"].keys():
        input.append(f"results/{dataset}/filtered/{dataset}.non_numts.tsv"),
        input.append(f"results/{dataset}/filtered/{dataset}.non_numts.clusters.fasta")
        taxonomy = config["input"][dataset]["taxonomy"]
        df = pd.read_csv(taxonomy, sep="\t", index_col=0)
        # Calculate the number of clusters for each order
        repdf = df.loc[df.representative==1]
        order_size = repdf.groupby("Order").size()
        # Only include orders with more than one cluster
        orders = order_size.loc[order_size>1].index.tolist()
        input += expand("results/{dataset}/evaluation/{order}_numt_{f}.tsv", dataset=[dataset], order=orders, f=["analysis", "evaluation"])
    return input

def read_orders(config):
    # generate tool-specific list of orders
    orders = {}

orders = read_orders(config)

rule numts:
    input:
        numts_input,

rule generate_order_seqs:
    """
    Generate subsets of the cluster representative ASV sequences for an order.
    """
    output:
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/numts_filtering/alignments/{order}_seqs.fasta"
    input:
        #taxonomy = lambda wildcards: config["input"][wildcards.dataset]["taxonomy"],
        taxonomy = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_taxonomy.tsv"
        #fasta = lambda wildcards: config["input"][wildcards.dataset]["fasta"],
        fasta = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_reps.fasta"
    log:
        #"results/{dataset}/alignments/generate_{order}_seqs.log"
        "logs/numts/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{order}_generate_order_seqs.log",
    conda: "../envs/r-env.yml"
    script: "../scripts/numt_filtering/generate_order_seqs.R"

rule generate_trimmed_seq:
    """
    Remove first base in ASV sequences
    """
    output:
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/numts_filtering/trimmed/{order}_seqs_trimmed.fasta"
    input:
        fasta=rules.generate_order_seqs.output,
    log:
        "logs/numts/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/{run_name}/{order}_generate_trimmed_seq.log",
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        indir=lambda wildcards, input: os.path.dirname(input.fasta),
    conda: "../envs/r-env.yml"
    script: "../scripts/numt_filtering/generate_trimmed_seqs.R"

rule generate_aa_seqs:
    """
    Translate to protein using code table 5 (invertebrate mt DNA)
    """
    output:
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/numts_filtering/aa/{order}_seqs_aa.fasta"
    input:
        fasta=rules.generate_trimmed_seq.output,
    log:
        "logs/numts/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{order}_generate_aa_seqs.log",
    conda: "../envs/r-env.yml"
    script: "../scripts/numt_filtering/generate_aa_seqs.R"

rule mafft_align:
    """
    Align protein sequences using MAFFT
    """
    output:
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/numts_filtering/mafft/{order}_seqs_aa_aligned.fasta"
    input:
        rules.generate_aa_seqs.output[0]
    log:
        "logs/numts/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{order}_mafft_align.log",
    conda: "../envs/mafft.yml"
    resources:
        runtime = 60 * 2,
        mem_mb = mem_allowed,
    threads: 4
    shell:
        "mafft --auto --thread {threads} {input} > {output} 2>{log}"

rule pal2nal:
    """
    Generate the corresponding nucleotide alignments with pal2nal
    """
    output:
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/numts_filtering/alignments/{order}_seqs_pal2nal.fasta"
    input:
        pep=rules.mafft_align.output[0],
        nuc=rules.generate_trimmed_seq.output[0],
    log:
        "logs/numts/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{order}_pal2nal.log"
    conda:
        "../envs/pal2nal.yml"
    params:
        table = "5"
    threads: 1
    resources:
        runtime = 60 * 2,
        mem_mb = mem_allowed,
    shell:
        """
        pal2nal.pl {input.pep} {input.nuc} -output fasta -codontable {params.table} > {output} 2>{log}
        """

rule combined_filter:
    """
    Filter the data for each order using with default settings.
    """
    output:
        tsv="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/numts_filtering/combined_filter/{order}_echo_analysis.tsv"
    input:
        fasta=rules.pal2nal.output[0],
        taxonomy = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_taxonomy.tsv"
        counts = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_counts.tsv"
    log:
        "logs/numts/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{order}_combined_filter.log"
    threads: lambda wildcards: 20 if wildcards.order in ["unclassified", "Diptera", "Hymenoptera", "unclassified.Arthropoda", "unclassified.Insecta"] else 1
    resources:
        runtime = 60*24, #lambda wildcards: 60 * 24 if wildcards.order in ["unclassified", "Diptera", "Hymenoptera"] else 60 * 4,
        mem_mb = mem_allowed,
    params:
        spikeins=config["numts"]["spikein_clusters"],
        functions="../workflow/scripts/numt_filtering/functions.R",
        codon_model="../workflow/scripts/numt_filtering/codon_model.R",
        n_closest=config["numts"]["n_closest"]
    conda: "../envs/r-env.yml"
    script: "../scripts/numt_filtering/combined_filter.R"    

rule abundance_filter:
    """
    Filter the data by abundance for an order.
    """
    output:
        tsv="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/numts_filtering/abundance_filter/{order}_abundance_analysis.tsv"
    input:
        taxonomy = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_taxonomy.tsv"
        counts = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_counts.tsv"
    log:
        "logs/numts/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{order}_abundance_filter.log"
    params:
        functions="../workflow/scripts/numt_filtering/functions.R",
    conda: "../envs/r-env.yml"
    script: "../scripts/numt_filtering/abundance_filter.R"

rule generate_cluster_analysis:
    output:
        tsv="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/numts_filtering/cluster_analysis/{order}_cluster_analysis.tsv"
    input:
        taxonomy = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_taxonomy.tsv"
        counts = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_counts.tsv"
        consensus_taxonomy = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_consensus_taxonomy.tsv"
    log:
        "logs/numts/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{order}_generate_cluster_analysis.log"
    params:
        functions = "../workflow/scripts/numt_filtering/functions.R",
    conda: "../envs/r-env.yml"
    script: "../scripts/numt_filtering/generate_cluster_analysis.R"

rule evaluate_order:
    output:
        numt_res="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/numts_filtering/evaluation/{order}_numt_analysis.tsv",
        numt_eval="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/numts_filtering/evaluation/{order}_numt_evaluation.tsv"        
    input:
        filter1=rules.combined_filter.output.tsv,
        filter2=rules.abundance_filter.output.tsv,
        clust_file=rules.generate_cluster_analysis.output.tsv,
        taxonomy="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_taxonomy.tsv"
    params:
        trusted=config["numts"]["non_numt_ASVs"],
        functions="../workflow/scripts/numt_filtering/functions.R",
    log:
        "logs/numts/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{order}_evaluate_order.log"
    conda: "../envs/r-env.yml"
    script: "../scripts/numt_filtering/evaluate_order.R"

def filtered_input(wc):
    """
    Generate the input for the filtered rule.
    """
    # first check if the cluster taxonomy file exists already
    taxonomy = f"results/{wc.tool}/{wc.rundir}/{wc.chimera_run}/{wc.chimdir}/{wc.rank}/runs/{wc.run_name}/cluster_taxonomy.tsv"
    if os.path.exists(taxonomy):
        taxfile = taxonomy
    # if not, use the taxonomy file from the input directory
    else:
        taxfile = f"{config[wc.rundir]["taxonomy"]}/asv_taxonomy.tsv"
    df = pd.read_csv(taxfile, sep="\t", index_col=0)
    df.rename(columns = lambda x: x.lower(), inplace=True)
    orders = df["order"].unique().tolist()
    return expand("results/{dataset}/evaluation/{order}_numt_analysis.tsv", dataset=[wildcards.dataset], order=orders)

rule filtered:
    """
    Combines the numt results for each order and outputs a taxonomy table of all non-numt ASVs 
    as well as a fasta file of the non-numt clusters.
    """
    output:
        tsv="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/non_numts_clusters.tsv",
        fasta="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/non_numts_clusters.fasta",
    input:
        numt_files=filtered_input,
        fasta = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_reps.fasta",
        taxonomy = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_taxonomy.tsv"
    log: 
        "logs/numts/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/filtered.log"
    conda: "../envs/r-env.yml"
    script: "../scripts/numt_filtering/filter.R"