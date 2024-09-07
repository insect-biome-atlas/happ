localrules:
    filter_seqs,

def get_input_fasta(wildcards):
    """
    Determine the input fasta file for the dataset
    """
    if wildcards.chimdir == "raw":
        f = expand(
            "data/{rundir}/asv_seqs.fasta",
            rundir=wildcards.rundir,
        )
    else:
        f = expand(
            "results/chimera/{rundir}/filtered/{chimera_run}/{chimdir}/nonchimeras.fasta",
            rundir=wildcards.rundir,
            chimera_run=wildcards.chimera_run,
            chimdir=wildcards.chimdir,
        )
    return f[0]

def get_input_taxa(wildcards):
    """
    Return the taxonomy file for the dataset
    """
    # check if the taxonomy source is a file
    if os.path.isfile(config["taxonomy_source"]):
        return taxonomy_source
    return expand(
        "results/{taxonomy_source}/{rundir}/{taxonomy_source}.tsv",
        taxonomy_source=config["taxonomy_source"],
        rundir=config["rundir"],
    )[0]

checkpoint filter_seqs:
    """
    Checkpoint for first round of filtering. Takes as input the chimera-filtered fasta
    if chimera filtering is activated. Ensures that all ASVs are present in both the counts file
    and the sequence file.
    """
    input:
        counts="data/{rundir}/asv_counts.tsv",
        fasta=get_input_fasta,
        tax=get_input_taxa,
    output:
        directory("results/common/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa")
    log:
        "logs/filter_seqs/{rundir}/{chimera_run}/{chimdir}/{rank}.filter.log",
    params:
        split_rank=config["split_rank"],
    script:
        "../scripts/common.py"

## VSEARCH ALIGNMENTS ##
rule vsearch_align:
    input:
        fasta="results/common/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/asv_seqs.fasta.gz",
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
        runtime=60 * 24,
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

def get_vsearch_files(wildcards):
    checkpoint_dir = checkpoints.filter_seqs.get(
        rundir=config["rundir"], 
        chimera_run=config["chimera"]["run_name"], 
        chimdir=config["chimdir"], 
        rank=config["split_rank"]
        ).output[0]
    files = expand("results/vsearch/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/asv_seqs.dist.gz",
        rundir=config["rundir"],
        chimera_run=config["chimera"]["run_name"],
        chimdir=config["chimdir"],
        rank=config["split_rank"],
        tax=glob_wildcards(os.path.join(checkpoint_dir, "{tax}", "asv_seqs.fasta.gz")).tax
        )
    return files

rule vsearch:
    """
    vsearch pseudo-target
    """
    input:
        get_vsearch_files,
        