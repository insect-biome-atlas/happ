localrules:
    clean_fasta,
    clean_counts,
    precision_recall_cleaned

rule postprocessing:
    """
    Target rule for postprocess cleaning of ASVs
    """
    input:
        expand("results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/{f}",
                tool=config["software"], rundir=config["rundir"], chimera_run=config["chimera"]["run_name"], 
                chimdir=config["chimdir"], rank=config["split_rank"], run_name=config["run_name"], 
                f=["non_numts.cleaned.tsv", "non_numts_clusters.cleaned.fasta", "precision_recall.non_numts.cleaned.txt"]),

rule clean_asvs:
    """
    Cleans ASVs from NUMTs filtered data based on presence in blanks
    """
    output:
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/non_numts.cleaned.tsv"
    input:
        tsv=rules.filtered.output.tsv,
        counts=f"data/{config['rundir']}/asv_counts.tsv",
        meta=config["metadata_file"],
    log:
        "logs/postprocess/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/clean_asvs.log"
    params:
        max_blank_occurrence = config["max_blank_occurrence"],
        split_col = config["split_col"],
        blank_val = config["blank_val"],
        sample_type_col = config["sample_type_col"],
        sample_id_col = config["sample_id_col"]
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
    Cleans the NUMTs filtered fasta file based on output above
    """
    output:
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/non_numts_clusters.cleaned.fasta"
    input:
        tsv=rules.clean_asvs.output[0],
        fasta=rules.filtered.output.fasta,
    log:
        "logs/postprocess/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/clean_fasta.log"
    script: "../scripts/common.py"

rule clean_counts:
    """
    Cleans the cluster counts based on output above
    """
    output:
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/non_numts_clusters_counts.cleaned.tsv"
    input:
        tsv=rules.clean_asvs.output[0],
        counts="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_counts.tsv",
    log:
        "logs/postprocess/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/clean_counts.log"
    script: "../scripts/common.py"

rule precision_recall_cleaned:
    input:
        clust_file=rules.clean_asvs.output[0]
    output:
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/precision_recall.non_numts.cleaned.txt",
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/precision_recall.non_numts.cleaned.order.txt",
    log:
        "logs/postprocess/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/precision_recall.log",
    params:
        src="workflow/scripts/evaluate_clusters.py",
        eval_rank=config["evaluation_rank"],
    shell:
        """
        python {params.src} {input[0]} {input[0]} --rank {params.eval_rank} --order_level {output[1]} > {output[0]} 2>{log}
        """