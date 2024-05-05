import pandas as pd
import os


localrules:
    add_sums,
    nonchimeric_taxa,
    nonchimeric_orders,
    append_size,
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
wildcard_constraints:
    sample=f"({'|'.join(samples)})",
    
## CHIMERA DETECTION ##
#######################

## Utilities ##

rule chimera_filtering:
    input:
        expand(
            "results/chimera/{rundir}/filtered/{chimera_run}/{method}.{algo}/{f}",
            rundir=config["rundir"],
            chimera_run=config["chimera"]["run_name"],
            method=config["chimera"]["method"],
            algo=config["chimera"]["algorithm"],
            f=["nonchimeras.fasta", "chimeras.fasta", "orders.txt", f"{config['split_rank']}.txt"],
        ),
        expand(
            "results/settings/{rundir}/{chimera_run}/{chimdir}/{run_name}.{suff}",
            rundir=config["rundir"],
            chimera_run=config["chimera"]["run_name"],
            chimdir=config["chimdir"],
            run_name=config["run_name"],
            suff=["json", "cmd"],
        ),

def filter_nonchimeric_taxa(taxfile, nonchimeras, rank, output):
    """
    Filter the taxonomy file to remove taxa in which all sequences have been marked as chimeric
    """
    df = pd.read_csv(taxfile, sep="\t", index_col=0)
    nonchim = []
    with open(nonchimeras, 'r') as fhin:
        for line in fhin:
            line = line.rstrip()
            if line.startswith(">"):
                asv = line.lstrip(">")
                nonchim.append(asv)
    df = df.loc[nonchim]
    taxa = list(df[rank].unique())
    with open(output, 'w') as fhout:
        for t in taxa:
            fhout.write(f"{t}\n")

rule nonchimeric_taxa:
    """
    Filter the taxonomy file to remove taxa in which all sequences have been marked as chimeric
    """
    output:
        expand("results/chimera/{rundir}/filtered/{chimera_run}/{method}.{algo}/{split_rank}.txt",
            rundir=config["rundir"],
            chimera_run=config["chimera"]["run_name"],
            method=config["chimera"]["method"],
            algo=config["chimera"]["algorithm"],
            split_rank=config["split_rank"],
        ),
    input:
        taxfile=expand("data/{rundir}/asv_taxa.tsv", rundir=config["rundir"]),
        nonchimeras=expand("results/chimera/{rundir}/filtered/{chimera_run}/{method}.{algo}/nonchimeras.fasta",
            rundir=config["rundir"],
            chimera_run=config["chimera"]["run_name"],
            method=config["chimera"]["method"],
            algo=config["chimera"]["algorithm"],
        )
    params:
        split_rank=config["split_rank"]
    run:
        filter_nonchimeric_taxa(taxfile=input.taxfile[0], nonchimeras=input.nonchimeras[0], rank=params.split_rank, output=output[0])

rule nonchimeric_orders:
    """
    Filter the taxonomy file to remove orders in which all sequences have been marked as chimeric
    """
    output:
        expand("results/chimera/{rundir}/filtered/{chimera_run}/{method}.{algo}/orders.txt",
            rundir=config["rundir"],
            chimera_run=config["chimera"]["run_name"],
            method=config["chimera"]["method"],
            algo=config["chimera"]["algorithm"],
        ),
    input:
        taxfile=expand("data/{rundir}/asv_taxa.tsv", rundir=config["rundir"]),
        nonchimeras=expand("results/chimera/{rundir}/filtered/{chimera_run}/{method}.{algo}/nonchimeras.fasta",
            rundir=config["rundir"],
            chimera_run=config["chimera"]["run_name"],
            method=config["chimera"]["method"],
            algo=config["chimera"]["algorithm"],
        )
    params:
        ranks = config["ranks"]
    run:
        ranks = params.ranks
        # case-insensitive search for 'order' in ranks
        rank = ranks[[x.lower() for x in ranks].index("order")]
        filter_nonchimeric_taxa(taxfile=input.taxfile[0], nonchimeras=input.nonchimeras[0], rank=rank, output=output[0])

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
    threads: 10
    params:
        src="workflow/scripts/sum_counts.py",
        tmpdir="$TMPDIR/{rundir}.sum_counts",
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
        src="workflow/scripts/add_size_to_fastaheader.py",
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
        fasta=temp("data/{rundir}/samplewise/{sample}.fasta"),
    input:
        sums="data/{rundir}/samplewise/{sample}.sum.tsv",
        fasta="data/{rundir}/asv_seqs.fasta",
    log:
        "logs/chimeras/{rundir}/{sample}.add-sums.log",
    params:
        src="workflow/scripts/add_size_to_fastaheader.py",
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
        temp(expand(
            "data/{{rundir}}/samplewise/{sample}.sum.tsv",
            sample=samples,
        )),
    input:
        counts="data/{rundir}/asv_counts.tsv",
    log:
        "logs/chimeras/{rundir}/split_counts.log",
    params:
        src="workflow/scripts/split_counts.py",
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        tmpdir="$TMPDIR/{rundir}.split",
    shell:
        """
        exec &> {log}
        mkdir -p {params.tmpdir}
        cp {input.counts} {params.tmpdir}/counts.tsv
        python {params.src} {params.tmpdir}/counts.tsv {params.tmpdir} >{log} 2>&1
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
        src="workflow/scripts/filter_chimeras.py",
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
        src="workflow/scripts/filter_chimeras.py",
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
