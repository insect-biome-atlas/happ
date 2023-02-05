import os.path

localrules:
    filter_seqs,


rule filter_seqs:
    input:
        counts=expand("data/{rundir}/asv_counts.tsv", rundir=config["rundir"]),
        fasta=expand("data/{rundir}/asv_seqs.fasta", rundir=config["rundir"]),
        tax=expand("data/{rundir}/asv_taxa.tsv", rundir=config["rundir"]),
    output:
        total_counts="results/common/{rundir}/{tax}/total_counts.tsv",
        counts="results/common/{rundir}/{tax}/asv_counts.tsv.gz",
        fasta="results/common/{rundir}/{tax}/asv_seqs.fasta.gz",
    log:
        "logs/filter_seqs/{rundir}/{tax}.filter.log"
    params:
        split_rank=config["split_rank"],
        tmpdir=os.path.expandvars("$TMPDIR/{rundir}_{tax}_filter_seqs"),
        total_counts=os.path.expandvars("$TMPDIR/{rundir}_{tax}_filter_seqs/total_counts.tsv"),
        counts=os.path.expandvars("$TMPDIR/{rundir}_{tax}_filter_seqs/asv_counts.tsv.gz"),
        fasta=os.path.expandvars("$TMPDIR/{rundir}_{tax}_filter_seqs/asv_seqs.fasta.gz"),
    script:
        "../scripts/common.py"


rule filter:
    input:
        expand(
            "results/common/{rundir}/{tax}/{f}",
            rundir=config["rundir"],
            tax=taxa,
            f=["total_counts.tsv", "asv_counts.tsv.gz", "asv_seqs.fasta.gz"],
        ),


rule vsearch_align:
    input:
        fasta="results/common/{rundir}/{tax}/asv_seqs.fasta.gz",
    output:
        dist="results/vsearch/{rundir}/{tax}/asv_seqs.dist.gz",
    log:
        "logs/vsearch/{rundir}/{tax}/vsearch_align.log",
    params:
        dist="$TMPDIR/vsearch/{rundir}/{tax}/asv_seqs.dist",
        fasta="$TMPDIR/vsearch/{rundir}/{tax}/asv_seqs.fasta",
        tmpdir="$TMPDIR/vsearch/{rundir}/{tax}",
        id=config["vsearch"]["id"],
        iddef=config["vsearch"]["iddef"],
        query_cov=config["vsearch"]["query_cov"],
    threads: config["vsearch"]["threads"]
    conda:
        "../envs/opticlust.yml"
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
    input:
        expand(
            "results/vsearch/{rundir}/{tax}/asv_seqs.dist.gz",
            rundir=config["rundir"],
            tax=taxa,
        ),
