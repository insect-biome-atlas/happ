localrules: 
    generate_cluster_analysis,
    generate_order_seqs,
    generate_trimmed_seq,
    generate_aa_seqs,
    abundance_filter,
    evaluate_order,
    precision_recall_numts,
    filtered,
    order_otutab,
    lulu_filtering,
    precision_recall_lulu,
    evaluate_order_lulu,
    filtered_lulu

## Target rules

rule numts_filtering:
    """
    Target rule for filtering non-numt ASV clusters from the data.
    """
    input:
        expand("results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/{f}",
                    tool=config["software"], rundir=config["rundir"], chimera_run=config["chimera"]["run_name"], 
                    chimdir=config["chimdir"], rank=config["split_rank"], run_name=config["run_name"], 
                    f=["non_numts.tsv", "non_numts_clusters.fasta", "precision_recall.non_numts.txt"]),

rule lulu_filtering:
    input:
        expand("results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/{f}",
                    tool=config["software"], rundir=config["rundir"], chimera_run=config["chimera"]["run_name"], 
                    chimdir=config["chimdir"], rank=config["split_rank"], run_name=config["run_name"], 
                    f=["non_numts.lulu.tsv", "non_numts_clusters.lulu.fasta", "precision_recall.lulu.txt"])        

## General utility rules
rule generate_order_seqs:
    """
    Generate subsets of the cluster representative ASV sequences for an order.
    """
    output:
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/numts_filtering/alignments/{order}_seqs.fasta"
    input:
        taxonomy = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_taxonomy.tsv",
        fasta = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_reps.fasta"
    log:
        "logs/numts/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{order}_generate_order_seqs.log",
    conda: "../envs/r-env.yml"
    script: "../scripts/numt_filtering/generate_order_seqs.R"

## Echo-algorithm based filtering rules

rule generate_trimmed_seq:
    """
    Remove first base in ASV sequences
    """
    output:
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/numts_filtering/trimmed/{order}_seqs_trimmed.fasta"
    input:
        fasta=rules.generate_order_seqs.output[0],
    log:
        "logs/numts/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{order}_generate_trimmed_seq.log",
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
    shell:
        """
        pal2nal.pl {input.pep} {input.nuc} -output fasta -codontable {params.table} > {output} 2>{log}
        """

# If a spike_in_file exists, attempt to identify the spike in clusters from the data
if os.path.exists(config["numts"]["spikein_file"]):
    rule id_spikeins:
        """
        Identifies spikein clusters in the data. These clusters are used to calibrate the counts during filtering.
        """
        output:
            "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/numts_filtering/spikeins.tsv"
        input:
            taxonomy_file = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_taxonomy.tsv",
            counts_file = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_counts.tsv",
            spikein_file = config["numts"]["spikein_file"]
        log:
            "logs/numts/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/id_spikeins.log"
        params:
            method = config["numts"]["spikein_method"],
        conda: "../envs/r-env.yml"
        script: "../scripts/numt_filtering/id_spikeins.R"
# otherwise, create an empty file
else:
    rule id_spikeins:
        """
        Creates empty spikein file
        """
        output:
            temp(touch("results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/numts_filtering/spikeins.tsv"))

rule combined_filter:
    """
    Filter the data for each order using with default settings.
    """
    output:
        tsv="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/numts_filtering/combined_filter/{order}_echo_analysis.tsv"
    input:
        fasta=rules.pal2nal.output[0],
        taxonomy = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_taxonomy.tsv",
        counts = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_counts.tsv",
        spikeins = rules.id_spikeins.output,
    log:
        "logs/numts/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{order}_combined_filter.log"
    threads: 1
    resources:
        runtime = lambda wildcards: 60*24*7 if wildcards.order in config["numts"]["large_orders"] else 60*12,
        mem_mb = lambda wildcards: 100000 if wildcards.order in config["numts"]["large_orders"] else 20000,
        slurm_partition = lambda wildcards: config["long_partition"] if wildcards.order in config["numts"]["large_orders"] else config["default_partition"]
    params:
        functions="workflow/scripts/numt_filtering/functions.R",
        codon_model="workflow/scripts/numt_filtering/codon_model.R",
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
        taxonomy = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_taxonomy.tsv",
        counts = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_counts.tsv"
    log:
        "logs/numts/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{order}_abundance_filter.log"
    params:
        functions="workflow/scripts/numt_filtering/functions.R",
        abundance_threshold=config["numts"]["abundance_threshold"]
    conda: "../envs/r-env.yml"
    script: "../scripts/numt_filtering/abundance_filter.R"

rule generate_cluster_analysis:
    output:
        tsv="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/numts_filtering/cluster_analysis/{order}_cluster_analysis.tsv"
    input:
        taxonomy = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_taxonomy.tsv",
        counts = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_counts.tsv",
        consensus_taxonomy = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_consensus_taxonomy.tsv"
    log:
        "logs/numts/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{order}_generate_cluster_analysis.log"
    params:
        functions = "workflow/scripts/numt_filtering/functions.R",
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
        functions="workflow/scripts/numt_filtering/functions.R",
        lulu=False
    log:
        "logs/numts/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{order}_evaluate_order.log"
    conda: "../envs/r-env.yml"
    script: "../scripts/numt_filtering/evaluate_order.R"

def get_orders(wc):
    # first check if the cluster taxonomy file exists already
    taxonomy = f"results/{wc.tool}/{wc.rundir}/{wc.chimera_run}/{wc.chimdir}/{wc.rank}/runs/{wc.run_name}/cluster_taxonomy.tsv"
    if os.path.exists(taxonomy):
        taxfile = taxonomy
    # if not, use the taxonomy file from the input directory
    else:
        taxfile = f"data/{config['rundir']}/asv_taxa.tsv"
    df = pd.read_csv(taxfile, sep="\t", index_col=0)
    df.rename(columns = lambda x: x.lower(), inplace=True)
    orders = df["order"].unique().tolist()
    return orders

def filtered_input(wc):
    """
    Generate the input for the filtered rule.
    """
    orders = get_orders(wc)
    if config["numts"]["filter_unclassified_rank"].lower() == "order":
        orders = [order for order in orders if not order.startswith("unclassified")]
    return expand("results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/numts_filtering/evaluation/{order}_numt_analysis.tsv",
        tool = wc.tool, rundir=wc.rundir, chimera_run=wc.chimera_run, chimdir=wc.chimdir, rank=wc.rank, run_name=wc.run_name, order=orders)

rule filtered:
    """
    Combines the numt results for each order and outputs a taxonomy table of all non-numt ASVs 
    as well as a fasta file of the non-numt clusters.
    """
    output:
        tsv="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/non_numts.tsv",
        fasta="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/non_numts_clusters.fasta",
    input:
        numt_files=filtered_input,
        fasta = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_reps.fasta",
        taxonomy = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_taxonomy.tsv"
    log: 
        "logs/numts/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/filtered.log"
    params:
        filter_unclassified_rank = config["numts"]["filter_unclassified_rank"],
    conda: "../envs/r-env.yml"
    script: "../scripts/numt_filtering/filter.R"

rule precision_recall_numts:
    input:
        clust_file=rules.filtered.output.tsv
    output:
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/precision_recall.non_numts.txt",
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/precision_recall.non_numts.order.txt",
    log:
        "logs/numts/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/precision_recall.log",
    params:
        src="workflow/scripts/evaluate_clusters.py",
        eval_rank=config["evaluation_rank"],
    shell:
        """
        python {params.src} {input[0]} {input[0]} --rank {params.eval_rank} --order_level {output[1]} > {output[0]} 2>{log}
        """

## LULU filtering rules

rule lulu_matchlist:
    """
    Create a matchlist for sequences in an order using vsearch
    """
    output:
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/numts_filtering/lulu/orders/{order}/matchlist.tsv"
    input:
        rules.generate_order_seqs.output[0],
    log:
        "logs/numts/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{order}_lulu_matchlist.log"
    conda:
        "../envs/vsearch.yml"
    params:
        tmpdir = "$TMPDIR/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/{order}/lulu",
    threads: 1
    shell:
        """
        mkdir -p {params.tmpdir}
        cat {input} | sed 's/>.\+ />/g' > {params.tmpdir}/cluster_reps.fasta
        vsearch --usearch_global {params.tmpdir}/cluster_reps.fasta --db {params.tmpdir}/cluster_reps.fasta --self --id .84 --iddef 1 \
            --userout {output} -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10 --threads {threads} > {log} 2>&1
        rm -rf {params.tmpdir}
        """

rule order_otutab:
    """
    Generate cluster counts table for an order
    """
    output:
        otutab="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/numts_filtering/counts/{order}_cluster_counts.tsv"
    input:
        counts="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_counts.tsv",
        tax="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_taxonomy.tsv"
    run:
        import pandas as pd
        order = wildcards.order
        tax = pd.read_csv(input.tax, sep="\t", index_col=0)
        tax.rename(columns = lambda x: x.lower(), inplace=True)
        clusters = pd.DataFrame(tax.loc[(tax["order"]==order)&(tax.representative==1), "cluster"])
        clusters = clusters.reset_index().set_index("cluster")
        rename_dict = clusters.to_dict()
        counts = pd.read_csv(input.counts, sep="\t", index_col=0)
        counts = counts.loc[clusters.index]
        counts.rename(index=rename_dict["ASV"], inplace=True)
        counts.to_csv(output.otutab, sep="\t")

rule lulu:
    output:
        curated_table="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/numts_filtering/lulu/orders/{order}/curated_table.tsv",
        otu_map="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/numts_filtering/lulu/orders/{order}/otu_map.tsv",
        log="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/numts_filtering/lulu/orders/{order}/log.txt",
    input:
        matchlist=rules.lulu_matchlist.output[0],
        otutab=rules.order_otutab.output.otutab,
    log:
        "logs/lulu_filtering/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/{order}_lulu_filtering.log"
    conda:
        "../envs/lulu.yaml"
    resources:
        runtime = lambda wildcards: 60*24*3 if wildcards.order in config["lulu"]["large_orders"] else 60*2,
        mem_mb = lambda wildcards: 100000 if wildcards.order in config["lulu"]["large_orders"] else 20000,
        slurm_partition = lambda wildcards: config["long_partition"] if wildcards.order in config["lulu"]["large_orders"] else config["default_partition"]
    shadow: "minimal"
    params:
        minimum_ratio_type = config["lulu"]["minimum_ratio_type"],
        minimum_ratio = config["lulu"]["minimum_ratio"], 
        minimum_match = config["lulu"]["minimum_match"], 
        minimum_relative_cooccurence = config["lulu"]["minimum_relative_cooccurence"],
    script:
        "../scripts/lulu.R"

rule evaluate_order_lulu:
    output:
        numt_res="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/numts_filtering/lulu/evaluation/{order}_analysis.tsv",
        numt_eval="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/numts_filtering/lulu/evaluation/{order}_evaluation.tsv"        
    input:
        otu_map=rules.lulu.output.otu_map,
        clust_file=rules.generate_cluster_analysis.output.tsv,
        taxonomy="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_taxonomy.tsv"
    params:
        trusted=config["lulu"]["non_numt_ASVs"],
        functions="workflow/scripts/numt_filtering/functions.R",
        lulu=True
    log:
        "logs/lulu_filtering/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{order}_evaluate_order.log"
    conda: "../envs/r-env.yml"
    script: "../scripts/numt_filtering/evaluate_order.R"

def filtered_input_lulu(wc):
    """
    Generate the input for the LULU filtered rule.
    """
    orders = get_orders(wc)
    if config["lulu"]["filter_unclassified_rank"].lower() == "order":
        orders = [order for order in orders if not order.startswith("unclassified")]
    return expand("results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/numts_filtering/lulu/evaluation/{order}_analysis.tsv",
        tool = wc.tool, rundir=wc.rundir, chimera_run=wc.chimera_run, chimdir=wc.chimdir, rank=wc.rank, run_name=wc.run_name, order=orders)

rule filtered_lulu:
    """
    Combines the LULU results for each order and outputs a taxonomy table of all non-numt ASVs 
    as well as a fasta file of the non-numt clusters.
    """
    output:
        tsv="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/non_numts.lulu.tsv",
        fasta="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/non_numts_clusters.lulu.fasta",
    input:
        numt_files=filtered_input_lulu,
        fasta = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_reps.fasta",
        taxonomy = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_taxonomy.tsv"
    log: 
        "logs/lulu_filtering/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/filtered.log",
    params:
        filter_unclassified_rank = config["lulu"]["filter_unclassified_rank"],
    conda: "../envs/r-env.yml"
    script: "../scripts/numt_filtering/filter.R"

rule precision_recall_lulu:
    input:
        rules.filtered_lulu.output.tsv
    output:
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/precision_recall.lulu.txt",
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/precision_recall.lulu.order.txt",
    log:
        "logs/lulu_filtering/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/precision_recall.log",
    params:
        src="workflow/scripts/evaluate_clusters.py",
        eval_rank=config["evaluation_rank"],
    shell:
        """
        python {params.src} {input[0]} {input[0]} --rank {params.eval_rank} --order_level {output[1]} > {output[0]} 2>{log}
        """