rule run_dbotu3:
    input:
        fasta = expand("data/{rundir}/asv_seqs.fasta", rundir=config["rundir"]),
        counts = expand("data/{rundir}/asv_counts.tsv", rundir=config["rundir"])
    output:
        "results/dbotu3/{rundir}/dbotu3.txt"
    log:
        log = "logs/dbotu3/{rundir}/dbotu3.log",
        err = "logs/dbotu3/{rundir}/dbotu3.err"
    params:
        settings = config["dbotu3"]["settings"]
    conda: "../envs/dbotu3.yml"
    shell:
        """
        dbotu3.py {input.counts} {input.fasta} -o {output} --log {log.log} 2> {log.err}
        """

rule dbotu3:
    input:
        expand("results/dbotu3/{rundir}/dbotu3.txt", rundir = config["rundir"])