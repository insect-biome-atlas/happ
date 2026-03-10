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
    taxonomy_source = config["taxonomy_source"]
    taxfiles = {
        "epa-ng": [
            config["epa-ng"]["msa"],
            config["epa-ng"]["tree"],
            config["epa-ng"]["ref_taxonomy"],
        ],
        "sintax": [config["sintax"]["ref"]],
        "vsearch": [config["qiime2"]["ref"], config["qiime2"]["taxfile"]],
    }
    # check if the taxonomy source is a link
    if os.path.islink(taxonomy_source):
        taxonomy_source_target = os.readlink(taxonomy_source)
        if os.path.isfile(taxonomy_source_target):
            return taxonomy_source
        else:
            raise FileNotFoundError(
                f"Taxonomy source {taxonomy_source_target} is not a file"
            )
    # check if the taxonomy source is a file
    elif os.path.isfile(taxonomy_source):
        return taxonomy_source
    elif taxonomy_source == "sintax+epa-ng" and all(
        os.path.exists(x) for x in taxfiles["epa-ng"] + taxfiles["sintax"]
    ):
        return expand(
            "results/taxonomy/sintax_epang/{rundir}/{heur}/taxonomy.tsv",
            rundir=config["rundir"],
            heur=config["epa-ng"]["heuristic"],
        )[0]
    elif taxonomy_source == "sintax" and all(
        os.path.exists(x) for x in taxfiles["sintax"]
    ):
        return expand(
            "results/taxonomy/{taxonomy_source}/{rundir}/taxonomy.tsv",
            taxonomy_source=taxonomy_source,
            rundir=config["rundir"],
        )[0]
    elif (
        taxonomy_source == "vsearch"
        and all(os.path.exists(x) for x in taxfiles["vsearch"])
        or all(os.path.exists(x) for x in taxfiles["sintax"])
    ):
        return expand(
            "results/taxonomy/vsearch/{rundir}/taxonomy.tsv",
            rundir=config["rundir"],
        )[0]
    elif taxonomy_source == "epa-ng" and all(
        os.path.exists(x) for x in taxfiles["epa-ng"]
    ):
        return expand(
            "results/taxonomy/epa-ng/{rundir}/{heur}/taxonomy.tsv",
            rundir=config["rundir"],
            heur=config["epa-ng"]["heuristic"],
        )[0]
    sys.exit(f"ERROR: Taxonomy reference files not found")


checkpoint split_input:
    """
    Splits the input fasta file into chunks
    """
    message:
        "Splitting input fasta for {wildcards.rundir}"
    output:
        directory("results/common/{rundir}/splits"),
    input:
        unpack(get_preprocessed_files),
    log:
        "logs/split_input/{rundir}.log",
    params:
        outdir=lambda wildcards, output: output[0],
        size=config["split_size"],
    resources:
        runtime=60,
    threads: 1
    shell:
        """
        cat {input.fasta} | seqkit split2 -O {params.outdir} -j {threads} -s {params.size} >{log} 2>&1
        """


checkpoint filter_seqs:
    """
    Checkpoint for first round of filtering. Takes as input the chimera-filtered fasta
    if chimera filtering is activated. Ensures that all ASVs are present in both the counts file
    and the sequence file.
    """
    message:
        "Generating fasta files per taxa at rank {wildcards.rank}"
    input:
        fasta=get_input_fasta,
        tax=get_input_taxa,
    output:
        directory("results/common/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa"),
    log:
        "logs/filter_seqs/{rundir}/{chimera_run}/{chimdir}/{rank}.filter.log",
    params:
        split_rank=config["split_rank"],
        src=workflow.source_path("../scripts/filter_seqs.py"),
    shadow:
        "minimal"
    threads: 1
    shell:
        """
        python {params.src} -f {input.fasta} -t {input.tax} -o {output} --tmpdir temp -j {threads} -r {params.split_rank} > {log} 2>&1
        """


rule filter_counts:
    input:
        fasta="results/common/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/asv_seqs.fasta.gz",
        counts="data/{rundir}/asv_counts.tsv",
    output:
        tsv="results/common/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/asv_counts.tsv.gz",
        tsv_total="results/common/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/total_counts.tsv",
    log:
        "logs/filter_counts/{rundir}/{chimera_run}/{chimdir}/{rank}/{tax}.filter.log",
    params:
        src=workflow.source_path("../scripts/filter_counts.py"),
    threads: 4
    shadow:
        "minimal"
    shell:
        """
        export POLARS_MAX_THREADS={threads}
        python {params.src} -f {input.fasta} -c {input.counts} -o .
        gzip -c asv_counts.tsv > {output.tsv}
        mv total_counts.tsv {output.tsv_total}
        """


## VSEARCH ALIGNMENTS ##
rule vsearch_align:
    message:
        "Aligning sequences for {wildcards.tax}"
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
        config["vsearch-env"]
    container:
        "docker://quay.io/biocontainers/vsearch:2.29.1--h6a68c12_0"
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
        rank=config["split_rank"],
    ).output[0]
    files = expand(
        "results/vsearch/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/asv_seqs.dist.gz",
        rundir=config["rundir"],
        chimera_run=config["chimera"]["run_name"],
        chimdir=config["chimdir"],
        rank=config["split_rank"],
        tax=glob_wildcards(
            os.path.join(checkpoint_dir, "{tax}", "asv_seqs.fasta.gz")
        ).tax,
    )
    return files


rule vsearch:
    """
    vsearch pseudo-target
    """
    input:
        get_vsearch_files,
