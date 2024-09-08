localrules: 
    generate_cluster_analysis,
    generate_order_seqs,
    generate_trimmed_seq,
    generate_aa_seqs,
    abundance_filter,
    evaluate_order,
    precision_recall_noise_filtered,
    filtered,
    noise_filtering,
    cleaning,
    order_otutab,
    lulu_filtering,
    precision_recall_lulu,
    evaluate_order_lulu,
    filtered_lulu

## Target rules

rule noise_filtering:
    """
    Target rule for filtering noise from the data.
    """
    input:
        expand("results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/runs/{noise_run}/{f}",
                    tool=config["software"], rundir=config["rundir"], chimera_run=config["chimera"]["run_name"], 
                    noise_run=config["noise_filtering"]["run_name"],
                    chimdir=config["chimdir"], rank=config["split_rank"], run_name=config["run_name"], 
                    f=["noise_filtered_cluster_taxonomy.tsv", "noise_filtered_clusters.fasta", "precision_recall.txt"]),

rule cleaning:
    """
    Target rule for cleaning of ASV clusters
    """
    input:
        expand("results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/runs/{noise_run}/cleaned/{f}",
                    tool=config["software"], rundir=config["rundir"], chimera_run=config["chimera"]["run_name"], 
                    noise_run=config["noise_filtering"]["run_name"],
                    chimdir=config["chimdir"], rank=config["split_rank"], run_name=config["run_name"], 
                    f=["cleaned_noise_filtered_cluster_taxonomy.tsv", "cleaned_noise_filtered_clusters.fasta", "precision_recall.txt"]),

rule lulu_filtering:
    input:
        expand("results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/lulu/runs/{lulu_run}/{f}",
                    tool=config["software"], rundir=config["rundir"], chimera_run=config["chimera"]["run_name"], 
                    chimdir=config["chimdir"], rank=config["split_rank"], run_name=config["run_name"],
                    lulu_run=config["lulu"]["run_name"],
                    f=["lulu_clusters.tsv", "lulu_clusters.fasta", "precision_recall.lulu.txt"])        


checkpoint split_clusters_by_taxa:
    """
    Checkpoint rule to split clusters by taxonomy prior to noise filtering.
    """
    output:
        directory("results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/cluster_analysis/")
    input:
        taxonomy = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_taxonomy.tsv",
        counts = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_counts.tsv",
        consensus_taxonomy = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_consensus_taxonomy.tsv",
        fasta = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_reps.fasta"
    log:
        "logs/noise_filtering/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/split_clusters_by_taxa.log"
    params:
        functions = "workflow/scripts/noise_filtering/functions.R",
    conda: "../envs/r-env.yml"
    script: "../scripts/noise_filtering/split_clusters_by_taxa.R"


rule generate_trimmed_seq:
    """
    Remove first base in ASV sequences
    """
    output:
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/trimmed/{order}_seqs_trimmed.fasta"
    input:
        fasta="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/cluster_analysis/{order}_seqs.fasta"
    log:
        "logs/noise_filtering/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{order}_generate_trimmed_seq.log",
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        indir=lambda wildcards, input: os.path.dirname(input.fasta),
    conda: "../envs/r-env.yml"
    script: "../scripts/noise_filtering/generate_trimmed_seqs.R"

rule generate_aa_seqs:
    """
    Translate to protein using code table 5 (invertebrate mt DNA)
    """
    output:
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/aa/{order}_seqs_aa.fasta"
    input:
        fasta=rules.generate_trimmed_seq.output,
    log:
        "logs/noise_filtering/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{order}_generate_aa_seqs.log",
    conda: "../envs/r-env.yml"
    script: "../scripts/noise_filtering/generate_aa_seqs.R"

rule mafft_align:
    """
    Align protein sequences using MAFFT
    """
    output:
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/mafft/{order}_seqs_aa_aligned.fasta"
    input:
        rules.generate_aa_seqs.output[0]
    log:
        "logs/noise_filtering/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{order}_mafft_align.log",
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
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/alignments/{order}_seqs_pal2nal.fasta"
    input:
        pep=rules.mafft_align.output[0],
        nuc=rules.generate_trimmed_seq.output[0],
    log:
        "logs/noise_filtering/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{order}_pal2nal.log"
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
if os.path.exists(config["noise_filtering"]["spikein_file"]):
    rule id_spikeins:
        """
        Identifies spikein clusters in the data. These clusters are used to calibrate the counts during filtering.
        """
        output:
            "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/spikeins.tsv"
        input:
            taxonomy_file = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_taxonomy.tsv",
            counts_file = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_counts.tsv",
            spikein_file = config["noise_filtering"]["spikein_file"]
        log:
            "logs/noise_filtering/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/id_spikeins.log"
        params:
            method = config["noise_filtering"]["spikein_method"],
        conda: "../envs/r-env.yml"
        script: "../scripts/noise_filtering/id_spikeins.R"
# otherwise, create an empty file
else:
    rule id_spikeins:
        """
        Creates empty spikein file
        """
        output:
            temp(touch("results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/spikeins.tsv"))

rule combined_filter:
    """
    Filter the data for each order using with default settings.
    """
    output:
        tsv="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/runs/{noise_run}/combined_filter/{order}_echo_analysis.tsv"
    input:
        fasta=rules.pal2nal.output[0],
        taxonomy = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_taxonomy.tsv",
        counts = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_counts.tsv",
        spikeins = rules.id_spikeins.output,
    log:
        "logs/noise_filtering/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/runs/{noise_run}/{order}_combined_filter.log"
    threads: 1
    resources:
        runtime = lambda wildcards: 60*24*3 if wildcards.order in config["noise_filtering"]["large_orders"] else 60,
        mem_mb = lambda wildcards: 100000 if wildcards.order in config["noise_filtering"]["large_orders"] else 20000,
    params:
        functions="workflow/scripts/noise_filtering/functions.R",
        codon_model="workflow/scripts/noise_filtering/codon_model.R",
        n_closest=config["noise_filtering"]["n_closest"],
        threshold=config["noise_filtering"]["threshold"],
        max_singleton_reads=config["noise_filtering"]["max_singleton_reads"],
        max_singletons=config["noise_filtering"]["max_singletons"],
    conda: "../envs/r-env.yml"
    script: "../scripts/noise_filtering/combined_filter.R"    

rule abundance_filter:
    """
    Filter the data by abundance for an order.
    """
    output:
        tsv="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/runs/{noise_run}/abundance_filter/{order}_abundance_analysis.tsv"
    input:
        taxonomy = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_taxonomy.tsv",
        counts = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_counts.tsv"
    log:
        "logs/noise_filtering/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/runs/{noise_run}/{order}_abundance_filter.log"
    params:
        functions="workflow/scripts/noise_filtering/functions.R",
        abundance_threshold=config["noise_filtering"]["abundance_threshold"]
    conda: "../envs/r-env.yml"
    script: "../scripts/noise_filtering/abundance_filter.R"

rule evaluate_order:
    output:
        noise_res="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/runs/{noise_run}/evaluation/{order}_noise_analysis.tsv",
        noise_eval="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/runs/{noise_run}/evaluation/{order}_noise_evaluation.tsv"        
    input:
        filter1=rules.combined_filter.output.tsv,
        filter2=rules.abundance_filter.output.tsv,
        clust_file=rules.generate_cluster_analysis.output.tsv,
        taxonomy="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_taxonomy.tsv"
    params:
        trusted=config["noise_filtering"]["non_noise_ASVs"],
        functions="workflow/scripts/noise_filtering/functions.R",
        lulu=False
    log:
        "logs/noise_filtering/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{noise_run}/{order}_evaluate_order.log"
    conda: "../envs/r-env.yml"
    script: "../scripts/noise_filtering/evaluate_order.R"

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
    if config["noise_filtering"]["filter_unclassified_rank"].lower() == "order":
        orders = [order for order in orders if not order.startswith("unclassified")]
    return expand("results/{{tool}}/{{rundir}}/{{chimera_run}}/{{chimdir}}/{{rank}}/runs/{{run_name}}/noise_filtering/runs/{{noise_run}}/evaluation/{order}_noise_analysis.tsv",
        order=orders)

rule filtered:
    """
    Combines the noise filtering results for each order and outputs a taxonomy table of all non-noisy ASVs 
    as well as a fasta file of the non-noisy clusters.
    """
    output:
        tsv="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/runs/{noise_run}/noise_filtered_cluster_taxonomy.tsv",
        fasta="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/runs/{noise_run}/noise_filtered_clusters.fasta",
    input:
        noise_files=filtered_input,
        fasta = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_reps.fasta",
        taxonomy = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_taxonomy.tsv"
    log: 
        "logs/noise_filtering/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{noise_run}/filtered.log"
    params:
        filter_unclassified_rank = config["noise_filtering"]["filter_unclassified_rank"],
    conda: "../envs/r-env.yml"
    script: "../scripts/noise_filtering/filter.R"

rule precision_recall_noise_filtered:
    output:
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/runs/{noise_run}/precision_recall.txt",
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/runs/{noise_run}/precision_recall.order.txt",
    input:
        clust_file=rules.filtered.output.tsv
    log:
        "logs/noise_filtering/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{noise_run}/precision_recall.log",
    params:
        src="workflow/scripts/evaluate_clusters.py",
        eval_rank=config["evaluation_rank"],
    shell:
        """
        python {params.src} {input[0]} {input[0]} --rank {params.eval_rank} --order_level {output[1]} > {output[0]} 2>{log}
        """

## Cleaning
if os.path.exists(config["cleaning"]["metadata_file"]):
    rule clean_asvs:
        """
        Cleans ASVs noise filtered data based on presence in blanks
        """
        output:
            "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/runs/{noise_run}/cleaned/cleaned_noise_filtered_cluster_taxonomy.tsv"
        input:
            tsv=rules.filtered.output.tsv,
            counts=f"data/{config['rundir']}/asv_counts.tsv",
            meta=config["cleaning"]["metadata_file"],
        log:
            "logs/cleaning/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/noise_filtering/{noise_run}/clean_asvs.log"
        params:
            max_blank_occurrence = config["cleaning"]["max_blank_occurrence"],
            split_col = config["cleaning"]["split_col"],
            blank_val = config["cleaning"]["blank_val"],
            sample_type_col = config["cleaning"]["sample_type_col"],
            sample_id_col = config["cleaning"]["sample_id_col"]
        resources:
            runtime = 120
        shell:
            """
            python workflow/scripts/clean_asvs.py --clustfile {input.tsv} --countsfile {input.counts} --metadata {input.meta} \
                --max_blank_occurrence {params.max_blank_occurrence} --split_col {params.split_col} --blank_val {params.blank_val} \
                --sample_type_col {params.sample_type_col} --sample_id_col {params.sample_id_col}  --outfile {output} 2>{log}
            """

    rule clean_fasta:
        """
        Cleans the noise filtered fasta file based on output above
        """
        output:
            "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/runs/{noise_run}/cleaned/cleaned_noise_filtered_clusters.fasta"
        input:
            tsv=rules.clean_asvs.output[0],
            fasta=rules.filtered.output.fasta,
        log:
            "logs/cleaning/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/noise_filtering/{noise_run}/clean_fasta.log"
        script: "../scripts/common.py"

    rule clean_counts:
        """
        Cleans the cluster counts based on output above
        """
        output:
            "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/runs/{noise_run}/cleaned_noise_filtered_cluster_counts.tsv"
        input:
            tsv=rules.clean_asvs.output[0],
            counts="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_counts.tsv",
        log:
            "logs/cleaning/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/noise_filtering/{noise_run}/clean_counts.log"
        script: "../scripts/common.py"

    rule precision_recall_cleaned:
        input:
            clust_file=rules.clean_asvs.output[0]
        output:
            "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/runs/{noise_run}/cleaned/precision_recall.txt",
            "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/runs/{noise_run}/cleaned/precision_recall.order.txt",
        log:
            "logs/cleaning/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/{noise_run}/precision_recall.log",
        params:
            src="workflow/scripts/evaluate_clusters.py",
            eval_rank=config["evaluation_rank"],
        shell:
            """
            python {params.src} {input[0]} {input[0]} --rank {params.eval_rank} --order_level {output[1]} > {output[0]} 2>{log}
            """

## LULU filtering rules
if config["lulu"]["alignment_tool"] == "vsearch":
    rule lulu_matchlist_vsearch:
        """
        Create a matchlist for sequences in an order using vsearch
        """
        output:
            "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/lulu/runs/{lulu_run}/orders/{order}/matchlist.tsv"
        input:
            rules.generate_order_seqs.output[0],
        log:
            "logs/lulu/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{lulu_run}/{order}_lulu_matchlist.log"
        conda:
            "../envs/vsearch.yml"
        params:
            tmpdir = "$TMPDIR/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/{order}/lulu",
            maxhits = config["lulu"]["max_target_seqs"]
        threads: 1
        shell:
            """
            mkdir -p {params.tmpdir}
            cat {input} | sed 's/>.\+ />/g' > {params.tmpdir}/cluster_reps.fasta
            vsearch --usearch_global {params.tmpdir}/cluster_reps.fasta --db {params.tmpdir}/cluster_reps.fasta --self --id .84 --iddef 1 \
                --userout {output} -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits {params.maxhits} --threads {threads} > {log} 2>&1
            rm -rf {params.tmpdir}
            """
else:
    rule lulu_matchlist_blastn:
        """
        Create a matchlist for sequences in an order using blastn
        """
        output:
            "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/lulu/runs/{lulu_run}/orders/{order}/matchlist.tsv"
        input:
            rules.generate_order_seqs.output[0],
        log:
            "logs/lulu/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{lulu_run}/{order}_lulu_matchlist.log"
        conda:
            "../envs/blastn.yml"
        params:
            tmpdir = "$TMPDIR/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/{order}/lulu",
            max_target_seqs = config["lulu"]["max_target_seqs"]
        threads: 1
        resources:
            runtime = lambda wildcards: 60*5 if wildcards.order in config["lulu"]["large_orders"] else 60,
            mem_mb = lambda wildcards: 100000 if wildcards.order in config["lulu"]["large_orders"] else 20000,
        shell:
            """
            mkdir -p {params.tmpdir}
            cat {input} | sed 's/>.\+ />/g' > {params.tmpdir}/cluster_reps.fasta
            makeblastdb -in {params.tmpdir}/cluster_reps.fasta -parse_seqids -dbtype nucl
            blastn -db {params.tmpdir}/cluster_reps.fasta -outfmt '6 qseqid sseqid pident' -out {output[0]} \
                -qcov_hsp_perc 80 -perc_identity 84 -query {params.tmpdir}/cluster_reps.fasta \
                -num_threads {threads} -max_target_seqs {params.max_target_seqs} 2>{log}
            rm -rf {params.tmpdir}
            """

rule order_otutab:
    """
    Generate cluster counts table for an order
    """
    output:
        otutab="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/counts/{order}_cluster_counts.tsv"
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
        curated_table="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/lulu/runs/{lulu_run}/orders/{order}/curated_table.tsv",
        otu_map="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/lulu/runs/{lulu_run}/orders/{order}/otu_map.tsv",
        log="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/lulu/runs/{lulu_run}/orders/{order}/log.txt",
    input:
        matchlist="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/lulu/runs/{lulu_run}/orders/{order}/matchlist.tsv",
        otutab=rules.order_otutab.output.otutab,
    log:
        "logs/lulu/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{lulu_run}/{order}_lulu_filtering.log"
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
        noise_res="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/lulu/runs/{lulu_run}/evaluation/{order}_analysis.tsv",
        noise_eval="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/lulu/runs/{lulu_run}/evaluation/{order}_evaluation.tsv"        
    input:
        otu_map=rules.lulu.output.otu_map,
        clust_file=rules.generate_cluster_analysis.output.tsv,
        taxonomy="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_taxonomy.tsv"
    params:
        trusted=config["lulu"]["non_noise_ASVs"],
        functions="workflow/scripts/noise_filtering/functions.R",
        lulu=True
    log:
        "logs/lulu/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{lulu_run}/{order}_evaluate_order.log"
    conda: "../envs/r-env.yml"
    script: "../scripts/noise_filtering/evaluate_order.R"

def filtered_input_lulu(wc):
    """
    Generate the input for the LULU filtered rule.
    """
    orders = get_orders(wc)
    if config["lulu"]["filter_unclassified_rank"].lower() == "order":
        orders = [order for order in orders if not order.startswith("unclassified")]
    return expand("results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/lulu/runs/{lulu_run}/evaluation/{order}_analysis.tsv",
        tool = wc.tool, rundir=wc.rundir, chimera_run=wc.chimera_run, chimdir=wc.chimdir, rank=wc.rank, run_name=wc.run_name, lulu_run=wc.lulu_run, order=orders)

rule filtered_lulu:
    """
    Combines the LULU results for each order and outputs a taxonomy table of all filtered ASVs 
    as well as a fasta file of the filtered clusters.
    """
    output:
        tsv="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/lulu/runs/{lulu_run}/lulu_clusters.tsv",
        fasta="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/lulu/runs/{lulu_run}/lulu_clusters.fasta",
    input:
        noise_files=filtered_input_lulu,
        fasta = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_reps.fasta",
        taxonomy = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_taxonomy.tsv"
    log: 
        "logs/lulu/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{lulu_run}/filtered.log"
    params:
        filter_unclassified_rank = config["lulu"]["filter_unclassified_rank"],
    conda: "../envs/r-env.yml"
    script: "../scripts/noise_filtering/filter.R"

rule precision_recall_lulu:
    input:
        rules.filtered_lulu.output.tsv
    output:
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/lulu/runs/{lulu_run}/precision_recall.lulu.txt",
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/noise_filtering/lulu/runs/{lulu_run}/precision_recall.lulu.order.txt",
    log:
        "logs/lulu/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{lulu_run}/precision_recall.log"
    params:
        src="workflow/scripts/evaluate_clusters.py",
        eval_rank=config["evaluation_rank"],
    shell:
        """
        python {params.src} {input[0]} {input[0]} --rank {params.eval_rank} --order_level {output[1]} > {output[0]} 2>{log}
        """

