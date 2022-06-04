localrules:
    lulu

rule run_lulu:
    input:
        dist = "results/vsearch/{rundir}/asv_seqs.dist.gz",
        counts = "results/common/{rundir}/asv_counts.tsv.gz"
    output:
        curated_table = "results/lulu/{rundir}/otus.tsv",
        otu_map = "results/lulu/{rundir}/otu_map.tsv"
    log:
        progress = "logs/lulu/{rundir}/progress.txt",
        log = "logs/lulu/{rundir}/log.txt"
    shadow: "minimal"
    conda: "../envs/lulu.yml"
    params:
        dist = "$TMPDIR/lulu/{rundir}/asv_seqs.dist",
        tmpdir = "$TMPDIR/lulu/{rundir}",
        minimum_ratio_type = config["lulu"]["minimum_ratio_type"],
        minimum_ratio = config["lulu"]["minimum_ratio"],
        minimum_match = config["lulu"]["minimum_match"],
        minimum_relative_cooccurence = config["lulu"]["minimum_relative_cooccurence"]
    resources:
        runtime = 60 * 24 * 10
    script:
        "../scripts/lulu.R"

rule lulu:
    input:
        expand("results/lulu/{rundir}/{f}.tsv", rundir = config["rundir"], f = ["otus", "otu_map"])