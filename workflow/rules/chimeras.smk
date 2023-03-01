import pandas as pd
import os


localrules:
    chimera_samplewise,
    add_sums,
    append_size,
    split_counts,
    filter_samplewise_chimeras,

def fetch_samples(f):
    """
    Look at the header of the counts file to figure out what samples to split by
    """
    r = pd.read_csv(f, sep="\t", nrows=1, index_col=0)
    return [x.replace(" ", "_") for x in list(r.columns)]


samples = fetch_samples(f=f"data/{config['rundir']}/asv_counts.tsv")

## CHIMERA DETECTION ##
#######################

## Utilities ##

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


rule chimera_batchwise:
    input:
        fasta=rules.append_size.output.fasta,
    output:
        chim="results/chimera/{rundir}/{chimera_run}/batchwise.{algo}/chimeras.fasta",
        nochim="results/chimera/{rundir}/{chimera_run}/batchwise.{algo}/nonchimeras.fasta",
        border="results/chimera/{rundir}/{chimera_run}/batchwise.{algo}/borderline.fasta",
        uchimeout="results/chimera/{rundir}/{chimera_run}/batchwise.{algo}/uchimeout.txt",
    log:
        "logs/chimeras/{rundir}/{chimera_run}/batchwise/{algo}.log",
    conda:
        "../envs/vsearch.yml"
    threads: 1
    resources:
        runtime=60 * 24,
    params:
        algorithm="--{algo}",
        abskew=get_abskew,
        dn=config["chimera"]["dn"],
        mindiffs=config["chimera"]["mindiffs"],
        mindiv=config["chimera"]["mindiv"],
        minh=config["chimera"]["minh"],
    shell:
        """
        vsearch --dn {params.dn} --mindiffs {params.mindiffs} --mindiv {params.mindiv} --minh {params.minh} \
            {params.abskew} --chimeras {output.chim} --borderline {output.border} --nonchimeras {output.nochim} \
            {params.algorithm} {input.fasta} --uchimeout {output.uchimeout} >{log} 2>&1
        """

rule chimera_samplewise:
    input:
        expand(
            "results/chimera/{rundir}/{chimera_run}/samplewise.{algo}/samples/{sample}/{f}.gz",
            rundir=config["rundir"],
            chimera_run=config["chimera"]["run_name"],
            algo=config["chimera"]["algorithm"],
            sample=samples,
            f=[
                "chimeras.fasta",
                "nonchimeras.fasta",
                "borderline.fasta",
                "uchimeout.txt",
                "uchimealns.out",
            ],
        ),
        expand(
            "results/chimera/{rundir}/{chimera_run}/samplewise.{algo}/{f}",
            rundir=config["rundir"],
            chimera_run=config["chimera"]["run_name"],
            algo=config["chimera"]["algorithm"],
            f=["chimeras.tsv", "nonchimeras.fasta"],
        ),

rule split_all_counts:
    input:
        expand(
            "data/{rundir}/samplewise/{sample}.sum.tsv",
            rundir=config["rundir"],
            sample=samples,
        )

rule split_counts:
    output:
        expand(
            "data/{{rundir}}/samplewise/{sample}.sum.tsv",
            sample=samples,
        ),
    input:
        counts="data/{rundir}/asv_counts.tsv",
    log:
        "logs/chimeras/{rundir}/split_counts.log",
    params:
        src=srcdir("../scripts/split_counts.py"),
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        tmpdir="$TMPDIR/{rundir}.split",
    conda:
        "../envs/polars.yml"
    shell:
        """
        mkdir -p {params.tmpdir}
        cp {input.counts} {params.tmpdir}/counts.tsv
        python {params.src} {params.tmpdir}/counts.tsv {params.tmpdir} 2>{log}
        rm {params.tmpdir}/counts.tsv
        mv {params.tmpdir}/* {params.outdir}
        rm -rf {params.tmpdir}
        """


rule add_sums:
    output:
        fasta="data/{rundir}/samplewise/{sample}.fasta",
    input:
        sums="data/{rundir}/samplewise/{sample}.sum.tsv",
        fasta="data/{rundir}/asv_seqs.fasta",
    log:
        "logs/chimeras/{rundir}/{sample}.add-sums.log",
    params:
        src=srcdir("../scripts/add_size_to_fastaheader.py"),
        tmpdir="$TMPDIR/{rundir}.{sample}.addsums",
    shell:
        """
        mkdir -p {params.tmpdir}
        python {params.src} {input.fasta} {input.sums} > {params.tmpdir}/{wildcards.sample}.fasta 2>{log}
        mv {params.tmpdir}/{wildcards.sample}.fasta {output.fasta}
        rm -rf {params.tmpdir}
        """


rule sample_chimera:
    output:
        chim="results/chimera/{rundir}/{chimera_run}/samplewise.{algo}/samples/{sample}/chimeras.fasta.gz",
        nochim="results/chimera/{rundir}/{chimera_run}/samplewise.{algo}/samples/{sample}/nonchimeras.fasta.gz",
        border="results/chimera/{rundir}/{chimera_run}/samplewise.{algo}/samples/{sample}/borderline.fasta.gz",
        uchimeout="results/chimera/{rundir}/{chimera_run}/samplewise.{algo}/samples/{sample}/uchimeout.txt.gz",
        alns="results/chimera/{rundir}/{chimera_run}/samplewise.{algo}/samples/{sample}/uchimealns.out.gz",
    input:
        fasta=rules.add_sums.output.fasta,
    log:
        "logs/chimeras/{rundir}/{chimera_run}/samplewise.{algo}/samples/{sample}.log",
    conda:
        "../envs/vsearch.yml"
    threads: 4
    resources:
        runtime=60 * 4,
    params:
        tmpdir="$TMPDIR/{rundir}.{chimera_run}.{algo}.{sample}.chim",
        outdir=lambda wildcards, output: os.path.dirname(output.chim),
        algorithm="--{algo}",
        abskew=get_abskew,
        dn=config["chimera"]["dn"],
        mindiffs=config["chimera"]["mindiffs"],
        mindiv=config["chimera"]["mindiv"],
        minh=config["chimera"]["minh"],
    shell:
        """
        if [ ! -s {input.fasta} ]; then
            touch {output.nochim} {output.chim} {output.alns} {output.border} {output.uchimeout} {params.outdir}/NOSEQS
        else
            mkdir -p {params.tmpdir}
            vsearch --threads {threads} --dn {params.dn} --mindiffs {params.mindiffs} \
                --mindiv {params.mindiv} --minh {params.minh} {params.abskew} \
                --chimeras {params.tmpdir}/chimeras.fasta --borderline {params.tmpdir}/borderline.fasta \
                --nonchimeras {params.tmpdir}/nonchimeras.fasta --uchimealns {params.tmpdir}/uchimealns.out \
                --uchimeout {params.tmpdir}/uchimeout.txt {params.algorithm} {input.fasta}  >{log} 2>&1
            gzip {params.tmpdir}/*
            mv {params.tmpdir}/*.gz {params.outdir}
            rm -rf {params.tmpdir}
        fi
        """


rule filter_samplewise_chimeras:
    output:
        nonchims="results/chimera/{rundir}/{chimera_run}/filtered/{samplewise_method}.{algo}/nonchimeras.fasta",
        chimeras="results/chimera/{rundir}/{chimera_run}/filtered/{samplewise_method}.{algo}/chimeras.tsv"
    input:
        fasta="data/{rundir}/asv_seqs.fasta",
        chims=expand("results/chimera/{rundir}/{chimera_run}/samplewise.{algo}/samples/{sample}/chimeras.fasta.gz",
            rundir=config["rundir"],
            chimera_run=config["chimera"]["run_name"],
            algo=config["chimera"]["algorithm"],
            sample=samples,
        ),
    log:
        "logs/chimeras/{rundir}/{chimera_run}/filtered/{samplewise_method}.{algo}/filter_samplewise_chimeras.log"
    params:
        tmpdir="$TMPDIR/{rundir}.{chimera_run}.{samplewise_method}.{algo}.filterchims",
        src=srcdir("../scripts/filter_samplewise_chimeras.py"),
    shell:
        """
        mkdir -p {params.tmpdir}
        python {params.src} --chims {input.chims} --fasta {input.fasta} --chimeraids {output.chimeras} > {params.tmpdir}/nonchimeras.fasta 2>{log}
        mv {params.tmpdir}/nonchimeras.fasta {output.nonchims}
        rm -rf {params.tmpdir}
        """

rule filter_batchwise_chimeras:
    """
    The filter batchwise rule takes as input the results from batchwise chimera 
    detection using the algorithm specified in the config and outputs nonchimeras
    under either 'strict' or default criteria. 
    """
    output:
        nonchims="results/chimera/{rundir}/{chimera_run}/filtered/{batchwise_method}.{algo}/nonchimeras.fasta",
    input:
        fasta="data/{rundir}/asv_seqs.fasta",
        uchimeout="results/chimera/{rundir}/{chimera_run}/batchwise.{algo}/uchimeout.txt",
        counts="data/{rundir}/asv_counts.tsv"
    params:
        src=srcdir("../scripts/filter_chimeras.py"),
        method="{batchwise_method}"
    conda:
        "../envs/polars.yml"
    shell:
        """
        python {params.src} --method {params.method} {input.fasta} {input.uchimeout} {input.counts}
        """

