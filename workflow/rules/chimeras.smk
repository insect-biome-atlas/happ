import pandas as pd
import os


localrules:
    add_sums,
    append_size,
    filter_samplewise_chimeras,
    filter_batchwise_chimeras,
    filter_chimeras


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

def get_preprocessed_files(wildcards):
    rundir = config["rundir"]
    if config["preprocessing"]["filter_codons"]:
        fastafile = f"results/preprocess/{rundir}/ASV_codon_filtered.fna"
        countsfile = f"results/preprocess/{rundir}/ASV_codon_filtered.table.tsv"
    elif config["preprocessing"]["filter_length"]:
        fastafile = f"results/preprocess/{rundir}/ASV_length_filtered.fna"
        countsfile = f"results/preprocess/{rundir}/ASV_length_filtered.table.tsv"
    else:
        fastafile = f"data/{rundir}/asv_seqs.fasta"
        countsfile = f"data/{rundir}/asv_counts.tsv"
    return {"fasta": fastafile, "counts": countsfile}

# TODO: Add checks for samples with zero sum counts after preprocessing?
samples = fetch_samples(f=f"data/{config['rundir']}/asv_counts.tsv")
wildcard_constraints:
    sample=f"({'|'.join(samples)})",
    
## CHIMERA DETECTION ##
#######################

## Utilities ##

rule filter_chimeras:
    """
    Pseudo-rule to act as a target for filtering chimeras
    """
    input:
        expand(
            "results/chimera/{rundir}/filtered/{chimera_run}/{method}.{algo}/{f}",
            rundir=config["rundir"],
            chimera_run=config["chimera"]["run_name"],
            method=config["chimera"]["method"],
            algo=config["chimera"]["algorithm"],
            f=["nonchimeras.fasta", "chimeras.fasta"]
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
    message: "Summing counts for sequences"
    input:
        unpack(get_preprocessed_files)
    output:
        sums="results/common/{rundir}/asv_sum.tsv",
    log:
        "logs/sum_asvs/{rundir}.log",
    params:
        src=workflow.source_path("../scripts/sum_counts.py"),
        tmpdir="$TMPDIR/{rundir}.sum_counts",
    shell:
        """
        python {params.src} {input.counts} > {output.sums} 2>{log}
        """

rule append_size:
    """
    Adds size annotation for batchwise mode
    """
    message: "Adding sizes to sequence headers"
    input:
        unpack(get_preprocessed_files),
        sums=rules.sum_asvs.output.sums,
    output:
        fasta="results/common/{rundir}/asv_seqs_size.fasta",
    log:
        "logs/append_size/{rundir}.log",
    params:
        src=workflow.source_path("../scripts/add_size_to_fastaheader.py"),
        tmpdir="$TMPDIR/{rundir}.addsums",
    shell:
        """
        mkdir -p {params.tmpdir}
        cp {input.fasta} {params.tmpdir}/asv_seqs.fasta
        python {params.src} {params.tmpdir}/asv_seqs.fasta {input.sums} > {params.tmpdir}/asv_seqs_size.fasta 2>{log}
        mv {params.tmpdir}/asv_seqs_size.fasta {output.fasta}
        rm -rf {params.tmpdir}
        """

### BATCHWISE ###


rule chimera_batchwise:
    """
    Run uchime algorithm with vsearch on the full dataset directly
    The wildcard 'algo' can be 'uchime_denovo', 'uchime2_denovo' or 'uchime3_denovo'
    """
    message: "Running chimera detection using {wildcards.algo} algorithm"
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
    conda: config["vsearch-env"]
    container: "docker://quay.io/biocontainers/vsearch:2.29.1--h6a68c12_0"
    params:
        algorithm="--{algo}",
        abskew=get_abskew,
        dn=config["chimera"]["dn"],
        mindiffs=config["chimera"]["mindiffs"],
        mindiv=config["chimera"]["mindiv"],
        minh=config["chimera"]["minh"],
    threads: 4
    shell:
        """
        vsearch --threads {threads} --dn {params.dn} --mindiffs {params.mindiffs} --mindiv {params.mindiv} \
            --minh {params.minh} {params.abskew} --chimeras {output.chim} \
            --borderline {output.border} --nonchimeras {output.nochim} \
            --uchimeout {output.uchimeout} --uchimealns {output.uchimealns} \
            {params.algorithm} {input.fasta}  >{log} 2>&1
        """


### SAMPLEWISE ###

rule split_counts_samplewise:
    message: "Generating counts file for {wildcards.sample}"
    output:
        temp("results/common/{rundir}/samplewise/{sample}.sum.tsv"),
    input:
        unpack(get_preprocessed_files),
    log:
        "logs/chimeras/{rundir}/{sample}.split_counts.log"
    params:
        src=workflow.source_path("../scripts/split_counts_samplewise.py"),
        tmpdir="$TMPDIR/split.{rundir}.{sample}"
    shell:
        """
        mkdir -p {params.tmpdir}
        cp {input.counts} {params.tmpdir}/counts.tsv
        python {params.src} {params.tmpdir}/counts.tsv {wildcards.sample} > {output[0]} 2>{log}
        rm -rf {params.tmpdir}
        """

rule add_sums:
    """
    Adds size annotation to fasta headers for samplewise mode
    """
    message: "Adding sizes to fasta headers for {wildcards.sample}"
    output:
        fasta=temp("results/common/{rundir}/samplewise/{sample}.fasta"),
    input:
        unpack(get_preprocessed_files),
        sums=rules.split_counts_samplewise.output[0],
    log:
        "logs/chimeras/{rundir}/{sample}.add-sums.log",
    params:
        src=workflow.source_path("../scripts/add_size_to_fastaheader.py"),
        tmpdir="$TMPDIR/{rundir}.{sample}.addsums",
    shell:
        """
        mkdir -p {params.tmpdir}
        cp {input.fasta} {params.tmpdir}/asv_seqs.fasta
        python {params.src} {params.tmpdir}/asv_seqs.fasta {input.sums} > {params.tmpdir}/{wildcards.sample}.fasta 2>{log}
        mv {params.tmpdir}/{wildcards.sample}.fasta {output.fasta}
        rm -rf {params.tmpdir}
        """


rule chimera_samplewise:
    """
    Run chimera detection on each sample
    """
    message: "Running chimera detection on {wildcards.sample} using {wildcards.algo} algorithm"
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
    conda: config["vsearch-env"]
    container: "docker://quay.io/biocontainers/vsearch:2.29.1--h6a68c12_0"
    threads: 4
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

def get_min_frac_chimeric_samples(config):
    if "min_frac_chimeric_samples" in config["chimera"]:
        return f"--min_frac_chimeric_samples {config['chimera']['min_frac_chimeric_samples']}"
    else:
        return ""

rule filter_samplewise_chimeras:
    """
    Filter samplewise chimera results
    """
    message: "Filtering chimeras using samplewise method"
    output:
        nonchims="results/chimera/{rundir}/filtered/{chimera_run}/samplewise.{algo}/nonchimeras.fasta",
        chimeras="results/chimera/{rundir}/filtered/{chimera_run}/samplewise.{algo}/chimeras.fasta",
    input:
        unpack(get_preprocessed_files),
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
        src=workflow.source_path("../scripts/filter_chimeras.py"),
        min_chimeric_samples=config["chimera"]["min_chimeric_samples"],
        min_frac_chimeric_samples=get_min_frac_chimeric_samples(config),
        minh=config["chimera"]["minh"],
        mindiff=config["chimera"]["mindiffs"],
        mindiv=config["chimera"]["mindiv"],
        algorithm=config["chimera"]["algorithm"],
    shell:
        """
        mkdir -p {params.tmpdir}
        cp {input.fasta} {params.tmpdir}/asv_seqs.fasta
        python {params.src} --uchimeout {input.uchimeout} --fasta {params.tmpdir}/asv_seqs.fasta \
            --min_chimeric_samples {params.min_chimeric_samples} {params.min_frac_chimeric_samples} \
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
    message: "Filtering chimeras using batchwise method"
    output:
        nonchims="results/chimera/{rundir}/filtered/{chimera_run}/batchwise.{algo}/nonchimeras.fasta",
        chims="results/chimera/{rundir}/filtered/{chimera_run}/batchwise.{algo}/chimeras.fasta",
        uchimeout="results/chimera/{rundir}/filtered/{chimera_run}/batchwise.{algo}/uchimeout.tsv",
    input:
        unpack(get_preprocessed_files),
        uchimeout=rules.chimera_batchwise.output.uchimeout,
    log:
        "logs/chimeras/{rundir}/filtered/{chimera_run}/batchwise.{algo}/filter_batchwise_chimeras.log",
    params:
        tmpdir="$TMPDIR/{rundir}.{chimera_run}.{algo}.filterchims_batchwise",
        src=workflow.source_path("../scripts/filter_chimeras.py"),
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
