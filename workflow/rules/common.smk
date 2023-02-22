import os.path


localrules:
    filter_seqs,
    append_size,


def get_filter_input(wildcards):
    if wildcards.algo == "none":
        f = expand(
            "data/{rundir}/asv_seqs.fasta",
            rundir=wildcards.rundir,
        )
    else:
        if config["samplewise_chimeras"]:
            f = (
                expand(
                    "results/chimera/{rundir}/samplewise/{algo}/nonchimeras.fasta",
                    rundir=wildcards.rundir,
                    algo=config["chimera_algorithm"],
                )
            )
        else:
            f = (
                expand(
                    "results/chimera/{rundir}/{algo}/nonchimeras.fasta",
                    rundir=wildcards.rundir,
                    algo=config["chimera_algorithm"],
                ),
            )
    return f[0]


rule filter_seqs:
    input:
        counts=expand("data/{rundir}/asv_counts.tsv", rundir=config["rundir"]),
        fasta=get_filter_input,
        tax=expand("data/{rundir}/asv_taxa.tsv", rundir=config["rundir"]),
    output:
        total_counts="results/common/{rundir}/{algo}/{tax}/total_counts.tsv",
        counts="results/common/{rundir}/{algo}/{tax}/asv_counts.tsv.gz",
        fasta="results/common/{rundir}/{algo}/{tax}/asv_seqs.fasta.gz",
    log:
        "logs/filter_seqs/{rundir}/{algo}/{tax}.filter.log",
    params:
        split_rank=config["split_rank"],
        tmpdir=os.path.expandvars("$TMPDIR/{rundir}_{algo}_{tax}_filter_seqs"),
        total_counts=os.path.expandvars(
            "$TMPDIR/{rundir}_{algo}_{tax}_filter_seqs/total_counts.tsv"
        ),
        counts=os.path.expandvars(
            "$TMPDIR/{rundir}_{algo}_{tax}_filter_seqs/asv_counts.tsv.gz"
        ),
        fasta=os.path.expandvars(
            "$TMPDIR/{rundir}_{algo}_{tax}_filter_seqs/asv_seqs.fasta.gz"
        ),
    script:
        "../scripts/common.py"


rule filter:
    """
    Pseudo-target for the filtering part of the workflow
    """
    input:
        expand(
            "results/common/{rundir}/{algo}/{tax}/{f}",
            rundir=config["rundir"],
            algo=config["chimera_algorithm"],
            tax=taxa,
            f=["total_counts.tsv", "asv_counts.tsv.gz", "asv_seqs.fasta.gz"],
        ),


## CHIMERA DETECTION ##
rule sum_asvs:
    input:
        counts="data/{rundir}/asv_counts.tsv",
    output:
        sums="data/{rundir}/asv_sum.tsv",
    log:
        "logs/sum_asvs/{rundir}.log",
    resources:
        runtime=60,
    threads: 10
    params:
        src=srcdir("../scripts/sum_counts.py"),
    conda:
        "../envs/polars.yml"
    shell:
        """
        python {params.src} {input.counts} > {output.sums} 2>{log}
        """


rule append_size:
    input:
        sums=rules.sum_asvs.output.sums,
        fasta="data/{rundir}/asv_seqs.fasta",
    output:
        fasta="data/{rundir}/asv_seqs_size.fasta",
    log:
        "logs/append_size/{rundir}.log",
    params:
        src=srcdir("../scripts/add_size_to_fastaheader.py"),
    shell:
        """
        python {params.src} {input.fasta} {input.sums} > {output.fasta} 2>{log}
        """


def get_abskew(wildcards):
    try:
        abskew = config["vsearch"]["abskew"]
    except KeyError:
        if wildcards.algo in ["uchime2_denovo", "uchime_denovo"]:
            abskew = 2.0
        elif wildcards.algo == "uchime3_denovo":
            abskew = 16.0
    return f"--abskew {abskew}"


rule chimera_full:
    input:
        fasta=rules.append_size.output.fasta,
    output:
        chim="results/chimera/{rundir}/{algo}/chimeras.fasta",
        nochim="results/chimera/{rundir}/{algo}/nonchimeras.fasta",
        border="results/chimera/{rundir}/{algo}/borderline.fasta",
        uchimeout="results/chimera/{rundir}/{algo}/uchimeout.txt",
    log:
        "logs/chimeras/{rundir}.{algo}.log",
    conda:
        "../envs/vsearch.yml"
    threads: 1
    resources:
        runtime=60 * 24,
    params:
        algorithm="--{algo}",
        abskew=get_abskew,
        dn=config["vsearch"]["dn"],
        mindiffs=config["vsearch"]["mindiffs"],
        mindiv=config["vsearch"]["mindiv"],
        minh=config["vsearch"]["minh"],
    shell:
        """
        vsearch --dn {params.dn} --mindiffs {params.mindiffs} --mindiv {params.mindiv} --minh {params.minh} \
            {params.abskew} --chimeras {output.chim} --borderline {output.border} --nonchimeras {output.nochim} \
            {params.algorithm} {input.fasta} --uchimeout {output.uchimeout} >{log} 2>&1
        """


## VSEARCH ALIGNMENTS ##
rule vsearch_align:
    input:
        fasta="results/common/{rundir}/{algo}/{tax}/asv_seqs.fasta.gz",
    output:
        dist="results/vsearch/{rundir}/{algo}/{tax}/asv_seqs.dist.gz",
    log:
        "logs/vsearch/{rundir}/{algo}/{tax}/vsearch_align.log",
    params:
        dist="$TMPDIR/vsearch/{rundir}/{algo}/{tax}/asv_seqs.dist",
        fasta="$TMPDIR/vsearch/{rundir}/{algo}/{tax}/asv_seqs.fasta",
        tmpdir="$TMPDIR/vsearch/{rundir}/{algo}/{tax}",
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
    """
    vsearch pseudo-target
    """
    input:
        expand(
            "results/vsearch/{rundir}/{algo}/{tax}/asv_seqs.dist.gz",
            rundir=config["rundir"],
            algo=config["chimera_algorithm"],
            tax=taxa,
        ),
