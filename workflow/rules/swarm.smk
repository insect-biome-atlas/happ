localrules:
    format_swarm,
    swarm

rule format_swarm:
    input:
        fasta = "results/common/{rundir}/asv_seqs.fasta.gz",
        counts = "results/common/{rundir}/total_counts.tsv"
    output:
        fasta = "results/swarm/{rundir}/reformat.fasta.gz"
    script:
        "../scripts/swarm_utils.py"

rule run_swarm:
    """
    swarm only requires that ASVs abundances are appended to fasta headers and
    that only ASVs with total abundances >0 are included
    """
    input:
        rules.format_swarm.output.fasta
    output:
        "results/swarm/{rundir}/swarm.txt"
    log:
        "logs/swarm/{rundir}/swarm.log"
    params:
        settings = config["swarm"]["settings"],
        tmpdir = "$TMPDIR/swarm/{rundir}",
        fasta = "$TMPDIR/swarm/{rundir}/reformat.fasta",
        output = "$TMPDIR/swarm/{rundir}/swarm.txt"
    threads: config["swarm"]["threads"]
    conda: "../envs/swarm.yml"
    resources:
        runtime = 60 * 24 * 10
    shell:
        """
        mkdir -p {params.tmpdir}
        gunzip -c {input} > {params.fasta}
        swarm {params.settings} {params.fasta} -o {params.output} -t {threads} > {log} 2>&1
        mv {params.output} {output}
        rm -rf {params.tmpdir}
        """

rule swarm:
    input:
        expand("results/swarm/{rundir}/swarm.txt", rundir = config["rundir"])