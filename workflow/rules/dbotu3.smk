localrules:
    dbotu3

rule run_dbotu3:
    """
    dbotu3 requires that all sequences in the fasta are also in the counts table
    and vice versa. However there's no requirement that total counts for 
    sequences are > 0
    """
    input:
        fasta = "results/common/{rundir}/asv_seqs.fasta.gz",
        counts = "results/common/{rundir}/asv_counts.tsv.gz"
    output:
        "results/dbotu3/{rundir}/dbotu3.txt"
    log:
        log = "logs/dbotu3/{rundir}/dbotu3.log",
        err = "logs/dbotu3/{rundir}/dbotu3.err"
    params:
        settings = config["dbotu3"]["settings"],
        tmpdir = "$TMPDIR/dbotu3/{rundir}",
        fasta = "$TMPDIR/dbotu3/{rundir}/asv_seqs.fasta",
        counts = "$TMPDIR/dbotu3/{rundir}/asv_counts.tsv",
        out = "$TMPDIR/dbotu3/{rundir}/dbotu3.txt"
    conda: "../envs/dbotu3.yml"
    resources:
        runtime = 60 * 24 * 10
    shell:
        """
        mkdir -p {params.tmpdir}
        gunzip -c {input.fasta} > {params.fasta}
        gunzip -c {input.counts} > {params.counts}
        dbotu3.py {params.counts} {params.fasta} -o {params.out} --log {log.log} 2> {log.err}
        mv {params.out} {output}
        rm -rf {params.tmpdir}
        """

rule dbotu3:
    input:
        expand("results/dbotu3/{rundir}/dbotu3.txt", rundir = config["rundir"])