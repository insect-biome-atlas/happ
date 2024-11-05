localrules:
    filter_codons,
    filter_length,

rule filter_length:
    """
    Filter ASVs by length.
    """
    output:
        fasta=ensure("results/preprocess/{rundir}/ASV_length_filtered.fna", non_empty=True),
        counts=ensure("results/preprocess/{rundir}/ASV_length_filtered.table.tsv", non_empty=True)
    input:
        fasta="data/{rundir}/asv_seqs.fasta",
        counts="data/{rundir}/asv_counts.tsv"
    log:
        "logs/preprocess/{rundir}/filter_length.log"
    params:
        min_length=lambda wildcards: f"-m {config['preprocessing']['min_length']}" if config['preprocessing']["min_length"] else "",
        max_length=lambda wildcards: f"-M {config['preprocessing']['max_length']}" if config['preprocessing']["max_length"] else "",
        outdir=lambda wildcards, output: os.path.dirname(output["fasta"]),
        src=workflow.source_path("../scripts/filt_length.py")
    shadow: "minimal"
    shell:
        """
        python {params.src} -f {input.fasta} -t {input.counts} -p ASV_length {params.min_length} {params.max_length} > {log} 2>&1
        mv ASV_length_filtered.fna {output.fasta}
        mv ASV_length_filtered.table.tsv {output.counts}
        """

rule filter_codons:
    """
    Filter ASVs with in-frame stop codons.
    """
    output:
        fasta=ensure("results/preprocess/{rundir}/ASV_codon_filtered.fna", non_empty=True),
        counts=ensure("results/preprocess/{rundir}/ASV_codon_filtered.table.tsv", non_empty=True),
        l="results/preprocess/{rundir}/ASV_codon_filtered.list"
    input:
        fasta=lambda wildcards: rules.filter_length.output.fasta if config["preprocessing"]["filter_length"] else f"data/{wildcards.rundir}/asv_seqs.fasta",
        counts=lambda wildcards: rules.filter_length.output.counts if config["preprocessing"]["filter_length"] else f"data/{wildcards.rundir}/asv_counts.tsv"
    log:
        "logs/preprocess/{rundir}/filter_codons.log"
    params:
        stop_codons=config['preprocessing']["stop_codons"],
        start_position=config['preprocessing']["start_position"],
        end_position=lambda wildcards: f"-e {config['preprocessing']['end_position']}" if config['preprocessing']['end_position']>0 else "",
        outdir=lambda wildcards, output: os.path.dirname(output["fasta"]),
        src=workflow.source_path("../scripts/filt_codons.py")
    shadow: "minimal"
    shell:
        """
        python {params.src} -f {input.fasta} -t {input.counts} -p ASV_codon -x {params.stop_codons} -s {params.start_position} {params.end_position} > {log} 2>&1
        mv ASV_codon_filtered.fna {output.fasta}
        mv ASV_codon_filtered.table.tsv {output.counts}
        mv ASV_codon_filtered.list {output.l}
        """