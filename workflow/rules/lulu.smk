localrules:
    lulu,
    lulu2tab,


rule download_lulu:
    output:
        src="src/lulu/Functions.R",
    log:
        "logs/lulu/download_lulu.log",
    params:
        url="https://raw.githubusercontent.com/tobiasgf/lulu/master/R/Functions.R",
    shell:
        """
        curl -L -o {output.src} {params.url} > {log} 2>&1
        """


rule run_lulu:
    input:
        src=rules.download_lulu.output.src,
        dist=rules.vsearch_align.output.dist,
        counts=rules.filter_seqs.output.counts,
    output:
        curated_table="results/lulu/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/{run_name}/otus.tsv",
        otu_map="results/lulu/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/{run_name}/otu_map.tsv",
    log:
        progress="logs/lulu/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/{run_name}/progress.txt",
        log="logs/lulu/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/{run_name}/log.txt",
    shadow:
        "minimal"
    conda:
        "../envs/lulu.yml"
    params:
        dist=os.path.expandvars("$TMPDIR/lulu/{rundir}/{chimera_run}/{chimdir}/{rank}/{tax}/asv_seqs.dist"),
        tmpdir=os.path.expandvars("$TMPDIR/lulu/{rundir}/{chimera_run}/{chimdir}/{rank}/{tax}"),
        minimum_ratio_type=config["lulu"]["minimum_ratio_type"],
        minimum_ratio=config["lulu"]["minimum_ratio"],
        minimum_match=config["lulu"]["minimum_match"],
        minimum_relative_cooccurence=config["lulu"]["minimum_relative_cooccurence"],
    resources:
        runtime=60 * 24 * 10,
        mem_mb=mem_allowed,
    script:
        "../scripts/lulu.R"


rule lulu2tab:
    input:
        rules.run_lulu.output.otu_map,
    output:
        "results/lulu/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/{run_name}/asv_clusters.tsv",
    params:
        tmpdir=os.path.expandvars("$TMPDIR/lulu/{rundir}/{chimera_run}/{chimdir}/{rank}/{tax}/{run_name}"),
        out=os.path.expandvars("$TMPDIR/lulu/{rundir}/{chimera_run}/{chimdir}/{rank}/{tax}/{run_name}/asv_clusters.tsv"),
    script:
        "../scripts/lulu_utils.py"


rule lulu:
    input:
        expand(
            "results/lulu/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/{run_name}/{f}.tsv",
            rundir=config["rundir"],
            chimera_run=config["chimera"]["run_name"],
            chimdir=config["chimdir"],
            rank=config["split_rank"],
            run_name=config["run_name"],
            f=["asv_clusters"],
            tax=taxa,
        ),
