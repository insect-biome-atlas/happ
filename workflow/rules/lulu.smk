localrules:
    lulu,
    lulu2tab

rule download_lulu:
    output:
        src = "src/lulu/Functions.R"
    log:
        "logs/lulu/download_lulu.log"
    params:
        url = "https://raw.githubusercontent.com/tobiasgf/lulu/master/R/Functions.R"
    shell:
        """
        curl -L -o {output.src} {params.url} > {log} 2>&1
        """

rule run_lulu:
    input:
        src = rules.download_lulu.output.src,
        dist = "results/vsearch/{rundir}/{tax}/asv_seqs.dist.gz",
        counts = "results/common/{rundir}/{tax}/asv_counts.tsv.gz"
    output:
        curated_table = "results/lulu/{rundir}/{tax}/otus.tsv",
        otu_map = "results/lulu/{rundir}/{tax}/otu_map.tsv"
    log:
        progress = "logs/lulu/{rundir}/{tax}/progress.txt",
        log = "logs/lulu/{rundir}/{tax}/log.txt"
    shadow: "minimal"
    conda: "../envs/lulu.yml"
    params:
        dist = "$TMPDIR/lulu/{rundir}/{tax}/asv_seqs.dist",
        tmpdir = "$TMPDIR/lulu/{rundir}/{tax}",
        minimum_ratio_type = config["lulu"]["minimum_ratio_type"],
        minimum_ratio = config["lulu"]["minimum_ratio"],
        minimum_match = config["lulu"]["minimum_match"],
        minimum_relative_cooccurence = config["lulu"]["minimum_relative_cooccurence"]
    resources:
        runtime = 60 * 24 * 10
    script:
        "../scripts/lulu.R"

rule lulu2tab:
    input:
        rules.run_lulu.output.otu_map
    output:
        "results/lulu/{rundir}/{tax}/asv_clusters.tsv"
    params:
        tmpdir = "$TMPDIR/lulu/{rundir}/{tax}",
        out = "$TMPDIR/lulu/{rundir}/{tax}/asv_clusters.tsv"
    script:
        "../scripts/lulu_utils.py"

rule lulu:
    input:
        expand("results/lulu/{rundir}/{tax}/{f}.tsv",
            rundir = config["rundir"], f = ["asv_clusters"], tax=taxa)
