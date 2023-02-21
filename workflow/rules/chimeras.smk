import pandas as pd
import os


localrules:
    chimera_samplewise,
    add_sums,
    split_counts,
    filter_samplewise_chimeras,

def fetch_samples(f):
    """
    Look at the header of the counts file to figure out what samples to split by
    """
    r = pd.read_csv(f, sep="\t", nrows=1, index_col=0)
    return [x.replace(" ", "_") for x in list(r.columns)]


samples = fetch_samples(f=f"data/{config['rundir']}/asv_counts.tsv")

# Read fasta and counts file and split by sample


rule chimera_samplewise:
    input:
        expand(
            "results/chimera/{rundir}/samplewise/{algo}/{sample}/{f}.gz",
            rundir=config["rundir"],
            algo=config["chimera_algorithm"],
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
            "results/chimera/{rundir}/samplewise/{algo}/{f}",
            rundir=config["rundir"],
            algo=config["chimera_algorithm"],
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
        chim="results/chimera/{rundir}/samplewise/{algo}/{sample}/chimeras.fasta.gz",
        nochim="results/chimera/{rundir}/samplewise/{algo}/{sample}/nonchimeras.fasta.gz",
        border="results/chimera/{rundir}/samplewise/{algo}/{sample}/borderline.fasta.gz",
        uchimeout="results/chimera/{rundir}/samplewise/{algo}/{sample}/uchimeout.txt.gz",
        alns="results/chimera/{rundir}/samplewise/{algo}/{sample}/uchimealns.out.gz",
    input:
        fasta=rules.add_sums.output.fasta,
    log:
        "logs/chimeras/{rundir}/samplewise/{algo}/{sample}.log",
    conda:
        "../envs/vsearch.yml"
    threads: 4
    resources:
        runtime=60 * 4,
    params:
        tmpdir="$TMPDIR/{rundir}.{algo}.{sample}.chim",
        outdir=lambda wildcards, output: os.path.dirname(output.chim),
        algorithm="--{algo}",
        abskew=get_abskew,
        dn=config["vsearch"]["dn"],
        mindiffs=config["vsearch"]["mindiffs"],
        mindiv=config["vsearch"]["mindiv"],
        minh=config["vsearch"]["minh"],
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
        nonchims="results/chimera/{rundir}/samplewise/{algo}/nonchimeras.fasta",
        chimeras="results/chimera/{rundir}/samplewise/{algo}/chimeras.tsv"
    input:
        fasta="data/{rundir}/asv_seqs.fasta",
        chims=expand("results/chimera/{rundir}/samplewise/{algo}/{sample}/chimeras.fasta.gz",
            rundir=config["rundir"],
            algo=config["chimera_algorithm"],
            sample=samples,
        ),
    log:
        "logs/chimeras/{rundir}/samplewise/{algo}/filter_samplewise_chimeras.log"
    params:
        tmpdir="$TMPDIR/{rundir}.{algo}.filterchims",
        src=srcdir("../scripts/filter_samplewise_chimeras.py"),
    shell:
        """
        mkdir -p {params.tmpdir}
        python {params.src} --chims {input.chims} --fasta {input.fasta} --chimeraids {output.chimeras} > {params.tmpdir}/nonchimeras.fasta 2>{log}
        mv {params.tmpdir}/nonchimeras.fasta {output.nonchims}
        rm -rf {params.tmpdir}
        """



