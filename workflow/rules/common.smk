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
    taxonomy_source=config["taxonomy_source"]
    # check if the taxonomy source is a file
    if os.path.isfile(taxonomy_source):
        return taxonomy_source
    return expand(
        "results/taxonomy/{taxonomy_source}/{rundir}/taxonomy.tsv",
        taxonomy_source=taxonomy_source,
        rundir=config["rundir"],
    )[0]

checkpoint split_input:
    """
    Splits the input fasta file into chunks with at most 500 sequences
    """
    output:
        directory("results/common/{rundir}/splits")
    input:
        unpack(get_preprocessed_files)
    log:
        "logs/split_input/{rundir}.log"
    params:
        outdir=lambda wildcards, output: output[0],
        size=500,
    resources:
        runtime = 60,
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
    shadow: "minimal"
    run:
        import pandas as pd
        from Bio import SeqIO
        import polars as pl
        import logging
        import subprocess
        logging.basicConfig(
            filename=log[0],
            filemode="w",
            level=logging.INFO,
            format="%(asctime)s - %(message)s",
        )
        outdir=output[0]
        tmpdir=os.environ.get("TMPDIR", "/tmp")
        split_rank = params.split_rank
        # open counts file lazily
        counts = pl.scan_csv(input.counts, separator="\t")
        # get name of first column
        first_col = counts.columns[0]
        # read taxonomy file
        taxdf = pd.read_csv(input.tax, sep="\t", index_col=0)
        # index the fasta file
        records = SeqIO.index(input.fasta, "fasta")
        # slice the taxonomy dataframe to only include the ASVs present in the fasta file
        taxdf = taxdf.loc[records.keys()]
        # get unique taxa to split by
        taxa = list(taxdf[split_rank].unique())
        n_taxa = len(taxa)
        for i, tax in enumerate(taxa, start=1):
            logging.info(f"Processing {tax} ({i}/{n_taxa})")
            os.makedirs(f"{outdir}/{tax}", exist_ok=True)
            taxtmpdir = f"{tmpdir}/filter_seqs/{wildcards.rundir}/{wildcards.chimera_run}/{wildcards.chimdir}/{wildcards.rank}/{tax}"
            os.makedirs(taxtmpdir, exist_ok=True)
            seqs_out = f"{taxtmpdir}/asv_seqs.fasta"
            counts_out = f"{taxtmpdir}/asv_counts.tsv"
            total_counts_out = f"{taxtmpdir}/total_counts.tsv"
            # get asvs for this taxon
            asvs = list(taxdf[taxdf[split_rank] == tax].index)
            # write asvs to file
            with open(f"{tax}.ids", "w") as f:
                f.write("\n".join(asvs))
            logging.info(f"Found {len(asvs)} ASVs for {tax}")
            # get sequences matching asvs using seqkit
            cmd = ["seqkit", "grep", "-f",f"{tax}.ids","-n", "--quiet", input.fasta]
            seqkit_process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            with open(seqs_out, "wb") as fhout:
                pigz_process = subprocess.Popen(["pigz", "-c"], stdin=seqkit_process.stdout, stdout=fhout)
                output, error = pigz_process.communicate()
            os.rename(seqs_out, f"{outdir}/{tax}/asv_seqs.fasta.gz")
            # filter the counts file and write
            tax_counts = counts.filter(pl.col(first_col).is_in(asvs)).collect()
            tax_counts.write_csv(counts_out, separator="\t")
            subprocess.run(["pigz", counts_out])
            os.rename(f"{counts_out}.gz", f"{outdir}/{tax}/asv_counts.tsv.gz")
            # Write the total counts file
            total_counts = tax_counts.with_columns(total=pl.sum_horizontal(pl.col("*").exclude(first_col))).select(first_col, "total")
            total_counts.columns = ["Representative_Sequence", "total"]
            total_counts.write_csv(total_counts_out, separator="\t")
            os.rename(total_counts_out, f"{outdir}/{tax}/total_counts.tsv")


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
        