import os.path


localrules:
    filter_seqs,


def get_filter_input(wildcards):
    if wildcards.chimdir == "raw":
        f = expand(
            "data/{rundir}/asv_seqs.fasta",
            rundir=config["rundir"],
        )
    else:
        f = expand(
            "results/chimera/{rundir}/filtered/{chimera_run}/{chimdir}/nonchimeras.fasta",
            rundir=config["rundir"],
            chimera_run=config["chimera"]["run_name"],
            chimdir=config["chimdir"],
        )
    return f[0]


rule filter_seqs:
    input:
        counts=expand("data/{rundir}/asv_counts.tsv", rundir=config["rundir"]),
        fasta=get_filter_input,
        tax=expand("data/{rundir}/asv_taxa.tsv", rundir=config["rundir"]),
    output:
        total_counts="results/common/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/total_counts.tsv",
        counts="results/common/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/asv_counts.tsv.gz",
        fasta="results/common/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/asv_seqs.fasta.gz",
    log:
        "logs/filter_seqs/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}.filter.log",
    params:
        split_rank=config["split_rank"],
        tmpdir=os.path.expandvars(
            "$TMPDIR/{rundir}_{chimera_run}_{chimdir}_{rank}_{tax}_filter_seqs"
        ),
        total_counts=os.path.expandvars(
            "$TMPDIR/{rundir}_{chimera_run}_{chimdir}_{rank}_{tax}_filter_seqs/total_counts.tsv"
        ),
        counts=os.path.expandvars(
            "$TMPDIR/{rundir}_{chimera_run}_{chimdir}_{rank}_{tax}_filter_seqs/asv_counts.tsv.gz"
        ),
        fasta=os.path.expandvars(
            "$TMPDIR/{rundir}_{chimera_run}_{chimdir}_{rank}_{tax}_filter_seqs/asv_seqs.fasta.gz"
        ),
    script:
        "../scripts/common.py"


rule filter:
    """
    Pseudo-target for the filtering part of the workflow
    """
    input:
        expand(
            "results/common/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/{f}",
            rundir=config["rundir"],
            chimera_run=config["chimera"]["run_name"],
            chimdir=config["chimdir"],
            rank=config["split_rank"],
            tax=taxa,
            f=["total_counts.tsv", "asv_counts.tsv.gz", "asv_seqs.fasta.gz"],
        ),


## VSEARCH ALIGNMENTS ##
rule vsearch_align:
    input:
        fasta=rules.filter_seqs.output.fasta,
    output:
        dist="results/vsearch/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/asv_seqs.dist.gz",
    log:
        "logs/vsearch/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/vsearch_align.log",
    params:
        dist="$TMPDIR/vsearch/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/asv_seqs.dist",
        fasta="$TMPDIR/vsearch/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/asv_seqs.fasta",
        tmpdir="$TMPDIR/vsearch/{rundir}/{chimera_run}/{chimdir}//{rank}/taxa/{tax}",
        id=config["vsearch"]["id"],
        iddef=config["vsearch"]["iddef"],
        query_cov=config["vsearch"]["query_cov"],
    threads: config["vsearch"]["threads"]
    conda:
        "../envs/vsearch.yml"
    resources:
        runtime=60 * 24 * 10,
    shell:
        """
        mkdir -p {params.tmpdir}
        gunzip -c {input.fasta} > {params.fasta}
        vsearch --usearch_global {params.fasta} --db {params.fasta} --self \
            --userout {params.dist} -userfields query+target+id --maxaccepts 0 --maxrejects 0 \
            --id {params.id} --iddef {params.iddef}  --query_cov {params.query_cov} --threads {threads} > {log} 2>&1
        gzip {params.dist}
        mv {params.dist}.gz {output.dist} 
        """


rule vsearch:
    """
    vsearch pseudo-target
    """
    input:
        expand(
            "results/vsearch/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/asv_seqs.dist.gz",
            rundir=config["rundir"],
            chimdir=config["chimdir"],
            chimera_run=config["chimera"]["run_name"],
            rank=config["split_rank"],
            tax=taxa,
        ),
