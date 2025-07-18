localrules:
    dbotu3,
    dbotu32tab,


rule run_dbotu3:
    """
    dbotu3 requires that all sequences in the fasta are also in the counts table
    and vice versa. However there's no requirement that total counts for 
    sequences are > 0
    """
    message: "Running dbotu3 clustering on sequences in {wildcards.tax}"
    input:
        fasta="results/common/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/asv_seqs.fasta.gz",
        counts="results/common/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/asv_counts.tsv.gz",
    output:
        tsv="results/dbotu3/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/{run_name}/dbotu3.tsv",
        memb="results/dbotu3/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/{run_name}/dbotu3.clusters.tsv",
    log:
        log="logs/dbotu3/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/{run_name}/dbotu3.log",
        err="logs/dbotu3/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/{run_name}/dbotu3.err",
    params:
        dist=config["dbotu3"]["dist"],
        abund=config["dbotu3"]["abund"],
        pval=config["dbotu3"]["pval"],
        tmpdir="$TMPDIR/dbotu3/{rundir}/{chimera_run}/{chimdir}/{rank}/{tax}",
        fasta="$TMPDIR/dbotu3/{rundir}/{chimera_run}/{chimdir}/{rank}/{tax}/asv_seqs.fasta",
        counts="$TMPDIR/dbotu3/{rundir}/{chimera_run}/{chimdir}/{rank}/{tax}/asv_counts.tsv",
        tsv="$TMPDIR/dbotu3/{rundir}/{chimera_run}/{chimdir}/{rank}/{tax}/dbotu3.tsv",
        memb="$TMPDIR/dbotu3/{rundir}/{chimera_run}/{chimdir}/{rank}/{tax}/dbotu3.clusters.tsv",
    conda: config["dbotu3-env"]
    shell:
        """
        mkdir -p {params.tmpdir}
        gunzip -c {input.fasta} > {params.fasta}
        gunzip -c {input.counts} > {params.counts}
        dbotu3.py -d {params.dist} -a {params.abund} -p {params.pval} \
            {params.counts} {params.fasta} -o {params.tsv} --membership {params.memb} --log {log.log} 2> {log.err}
        mv {params.tsv} {output.tsv}
        mv {params.memb} {output.memb}
        rm -rf {params.tmpdir}
        """


rule dbotu32tab:
    message: "Generating dbotu3 cluster file for {wildcards.tax}"
    input:
        rules.run_dbotu3.output.memb,
    output:
        "results/dbotu3/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/{run_name}/asv_clusters.tsv",
    params:
        tmpdir=os.path.expandvars("$TMPDIR/dbotu3/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/{run_name}"),
        out=os.path.expandvars("$TMPDIR/dbotu3/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/{run_name}/asv_clusters.tsv"),
    script:
        "../scripts/dbotu3_utils.py"


def get_dbotu3_files(wildcards):
    checkpoint_dir = checkpoints.filter_seqs.get(
        rundir=config["rundir"], 
        chimera_run=config["chimera"]["run_name"], 
        chimdir=config["chimdir"], 
        rank=config["split_rank"]
        ).output[0]
    files = expand("results/dbotu3/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/{run_name}/asv_clusters.tsv",
        rundir=config["rundir"],
        chimera_run=config["chimera"]["run_name"],
        chimdir=config["chimdir"],
        rank=config["split_rank"],
        tax=glob_wildcards(os.path.join(checkpoint_dir, "{tax}", "asv_seqs.fasta.gz")).tax,
        run_name=config["run_name"]
        )
    return files

rule dbotu3:
    input:
        get_dbotu3_files,
