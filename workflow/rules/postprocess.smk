localrules:
    clean_fasta,
    clean_counts

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
