import pandas as pd
import os


localrules:
    chimera_samplewise,
    add_sums,
    append_size,
    split_counts,
    filter_samplewise_chimeras,
    filter_batchwise_chimeras,


def fetch_samples(f):
    """
    Look at the header of the counts file to figure out what samples to split by
    """
    r = pd.read_csv(f, sep="\t", nrows=1, index_col=0)
    return [x.replace(" ", "_") for x in list(r.columns)]


def get_abskew(wildcards):
    try:
        abskew = config["chimera"]["abskew"]
    except KeyError:
        if wildcards.algo in ["uchime2_denovo", "uchime_denovo"]:
            abskew = 2.0
        elif wildcards.algo == "uchime3_denovo":
            abskew = 16.0
    return f"--abskew {abskew}"


samples = fetch_samples(f=f"data/{config['rundir']}/asv_counts.tsv")

## CHIMERA DETECTION ##
#######################

## Utilities ##


rule chimera_filtering:
    input:
        expand(
            "results/chimera/{rundir}/filtered/{chimera_run}/{method}.{algo}/{f}.fasta",
            rundir=config["rundir"],
            chimera_run=config["chimera"]["run_name"],
            method=config["chimera"]["method"],
            algo=config["chimera"]["algorithm"],
            f=["nonchimeras", "chimeras"],
        ),
        expand(
            "results/settings/{rundir}/{chimera_run}/{chimdir}/{run_name}.{suff}",
            rundir=config["rundir"],
            chimera_run=config["chimera"]["run_name"],
            chimdir=config["chimdir"],
            run_name=config["run_name"],
            suff=["json", "cmd"],
        ),


rule sum_asvs:
    """
    Sums up counts for batchwise mode
    """
    input:
        counts="data/{rundir}/asv_counts.tsv",
    output:
        sums="data/{rundir}/asv_sum.tsv",
    log:
        "logs/sum_asvs/{rundir}.log",
    resources:
        runtime=60,
        mem_mb=mem_allowed,
    threads: 10
    params:
        src=srcdir("../scripts/sum_counts.py"),
        tmpdir="$TMPDIR/{rundir}.sum_counts",
    conda:
        "../envs/polars.yml"
    shell:
        """
        mkdir -p {params.tmpdir}
        cp {input.counts} {params.tmpdir}/asv_counts.tsv
        python {params.src} {params.tmpdir}/asv_counts.tsv > {params.tmpdir}/asv_sum.tsv 2>{log}
        mv {params.tmpdir}/asv_sum.tsv {output.sums}
        rm -rf {params.tmpdir}
        """


rule append_size:
    """
    Adds size annotation for batchwise mode
    """
    input:
        sums=rules.sum_asvs.output.sums,
        fasta="data/{rundir}/asv_seqs.fasta",
    output:
        fasta="data/{rundir}/asv_seqs_size.fasta",
    log:
        "logs/append_size/{rundir}.log",
    params:
        src=srcdir("../scripts/add_size_to_fastaheader.py"),
        tmpdir="$TMPDIR/{rundir}.addsums",
    shell:
        """
        mkdir -p {params.tmpdir}
        cp {input.fasta} {params.tmpdir}/asv_seqs.fasta
        python {params.src} {params.tmpdir}/asv_seqs.fasta {input.sums} > {params.tmpdir}/asv_seqs_size.fasta 2>{log}
        mv {params.tmpdir}/asv_seqs_size.fasta {output.fasta}
        rm -rf {params.tmpdir}
        """


rule add_sums:
    """
    Adds size annotation to fasta headers for samplewise mode
    """
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
        cp {input.fasta} {params.tmpdir}/asv_seqs.fasta
        python {params.src} {params.tmpdir}/asv_seqs.fasta {input.sums} > {params.tmpdir}/{wildcards.sample}.fasta 2>{log}
        mv {params.tmpdir}/{wildcards.sample}.fasta {output.fasta}
        rm -rf {params.tmpdir}
        """


### BATCHWISE ###


rule chimera_batchwise:
    """
    Run uchime algorithm with vsearch on the full dataset directly
    The wildcard 'algo' can be 'uchime_denovo', 'uchime2_denovo' or 'uchime3_denovo'
    """
    input:
        fasta=rules.append_size.output.fasta,
    output:
        chim="results/chimera/{rundir}/batchwise.{algo}/chimeras.fasta",
        nochim="results/chimera/{rundir}/batchwise.{algo}/nonchimeras.fasta",
        border="results/chimera/{rundir}/batchwise.{algo}/borderline.fasta",
        uchimeout="results/chimera/{rundir}/batchwise.{algo}/uchimeout.txt",
        uchimealns="results/chimera/{rundir}/batchwise.{algo}/uchimealns.out",
    log:
        "logs/chimeras/{rundir}/batchwise/{algo}.log",
    conda:
        "../envs/vsearch.yml"
    threads: 1
    resources:
        runtime=60 * 24,
        mem_mb=mem_allowed,
    params:
        algorithm="--{algo}",
        abskew=get_abskew,
        dn=config["chimera"]["dn"],
        mindiffs=config["chimera"]["mindiffs"],
        mindiv=config["chimera"]["mindiv"],
        minh=config["chimera"]["minh"],
    shell:
        """
        vsearch --dn {params.dn} --mindiffs {params.mindiffs} --mindiv {params.mindiv} \
            --minh {params.minh} {params.abskew} --chimeras {output.chim} \
            --borderline {output.border} --nonchimeras {output.nochim} \
            --uchimeout {output.uchimeout} --uchimealns {output.uchimealns} \
            {params.algorithm} {input.fasta}  >{log} 2>&1
        """


### SAMPLEWISE ###


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


rule chimera_samplewise:
    """
    Run chimera detection on each sample
    """
    output:
        chim="results/chimera/{rundir}/samplewise.{algo}/samples/{sample}/chimeras.fasta.gz",
        nochim="results/chimera/{rundir}/samplewise.{algo}/samples/{sample}/nonchimeras.fasta.gz",
        border="results/chimera/{rundir}/samplewise.{algo}/samples/{sample}/borderline.fasta.gz",
        uchimeout="results/chimera/{rundir}/samplewise.{algo}/samples/{sample}/uchimeout.txt.gz",
        alns="results/chimera/{rundir}/samplewise.{algo}/samples/{sample}/uchimealns.out.gz",
    input:
        fasta=rules.add_sums.output.fasta,
    log:
        "logs/chimeras/{rundir}/samplewise.{algo}/samples/{sample}.log",
    conda:
        "../envs/vsearch.yml"
    threads: 4
    resources:
        runtime=60 * 4,
        mem_mb=mem_allowed,
    params:
        tmpdir="$TMPDIR/{rundir}.{algo}.{sample}.chim",
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


### FILTERING ###


rule filter_samplewise_chimeras:
    """
    Filter samplewise chimera results
    """
    output:
        nonchims="results/chimera/{rundir}/filtered/{chimera_run}/samplewise.{algo}/nonchimeras.fasta",
        chimeras="results/chimera/{rundir}/filtered/{chimera_run}/samplewise.{algo}/chimeras.fasta",
    input:
        fasta="data/{rundir}/asv_seqs.fasta",
        uchimeout=expand(
            "results/chimera/{rundir}/samplewise.{algo}/samples/{sample}/uchimeout.txt.gz",
            rundir=config["rundir"],
            algo=config["chimera"]["algorithm"],
            sample=samples,
        ),
    log:
        "logs/chimeras/{rundir}/filtered/{chimera_run}/samplewise.{algo}/filter_samplewise_chimeras.log",
    params:
        tmpdir="$TMPDIR/{rundir}.{chimera_run}.{algo}.filterchims_samplewise",
        src=srcdir("../scripts/filter_chimeras.py"),
        min_chimeric_samples=config["chimera"]["min_chimeric_samples"],
        minh=config["chimera"]["minh"],
        mindiff=config["chimera"]["mindiffs"],
        mindiv=config["chimera"]["mindiv"],
        algorithm=config["chimera"]["algorithm"],
    shell:
        """
        mkdir -p {params.tmpdir}
        cp {input.fasta} {params.tmpdir}/asv_seqs.fasta
        python {params.src} --uchimeout {input.uchimeout} --fasta {params.tmpdir}/asv_seqs.fasta \
            --min_chimeric_samples {params.min_chimeric_samples} \
            --chimfasta {params.tmpdir}/chimeras.fasta --nonchimfasta {params.tmpdir}/nonchimeras.fasta \
            --mindiff {params.mindiff} --mindiv {params.mindiv} --algorithm {params.algorithm} \
            --minh {params.minh} 2>{log}
        mv {params.tmpdir}/chimeras.fasta {output.chimeras}
        mv {params.tmpdir}/nonchimeras.fasta {output.nonchims}
        rm -rf {params.tmpdir}
        """


rule filter_batchwise_chimeras:
    """
    The filter batchwise rule takes as input the results from batchwise chimera 
    detection using the algorithm specified in the config and outputs nonchimeras
    under criteria matching 'min_samples_shared' or 'min_frac_samples_shared'
    """
    output:
        nonchims="results/chimera/{rundir}/filtered/{chimera_run}/batchwise.{algo}/nonchimeras.fasta",
        chims="results/chimera/{rundir}/filtered/{chimera_run}/batchwise.{algo}/chimeras.fasta",
        uchimeout="results/chimera/{rundir}/filtered/{chimera_run}/batchwise.{algo}/uchimeout.tsv",
    input:
        fasta="data/{rundir}/asv_seqs.fasta",
        uchimeout=rules.chimera_batchwise.output.uchimeout,
        counts="data/{rundir}/asv_counts.tsv",
    log:
        "logs/chimeras/{rundir}/filtered/{chimera_run}/batchwise.{algo}/filter_batchwise_chimeras.log",
    params:
        tmpdir="$TMPDIR/{rundir}.{chimera_run}.{algo}.filterchims_batchwise",
        src=srcdir("../scripts/filter_chimeras.py"),
        min_samples_shared=config["chimera"]["min_samples_shared"],
        min_frac_samples_shared=config["chimera"]["min_frac_samples_shared"],
        minh=config["chimera"]["minh"],
        mindiff=config["chimera"]["mindiffs"],
        mindiv=config["chimera"]["mindiv"],
        algorithm=config["chimera"]["algorithm"]
    shell:
        """
        mkdir -p {params.tmpdir}
        cp {input.counts} {params.tmpdir}/asv_counts.tsv
        cp {input.fasta} {params.tmpdir}/asv_seqs.fasta
        python {params.src} --min_samples_shared {params.min_samples_shared} \
            --min_frac_samples_shared {params.min_frac_samples_shared} \
            --fasta {params.tmpdir}/asv_seqs.fasta --uchimeout {input.uchimeout} \
            --counts {params.tmpdir}/asv_counts.tsv --chimfasta {params.tmpdir}/chimeras.fasta \
            --nonchimfasta {params.tmpdir}/nonchimeras.fasta --algorithm {params.algorithm} \
            --filteredout {params.tmpdir}/uchimeout.tsv --minh {params.minh} \
            --mindiff {params.mindiff} --mindiv {params.mindiv} 2>{log}
        mv {params.tmpdir}/chimeras.fasta {output.chims}
        mv {params.tmpdir}/nonchimeras.fasta {output.nonchims}
        mv {params.tmpdir}/uchimeout.tsv {output.uchimeout}
        rm -rf {params.tmpdir}
        """
