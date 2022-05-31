rule format_swarm:
    input:
        fasta = expand("data/{rundir}/asv_seqs.fasta", rundir=config["rundir"]),
        counts = expand("data/{rundir}/asv_counts.tsv", rundir=config["rundir"])
    output:
        fasta = "results/swarm/{rundir}/reformat.fasta"
    script:
        "../scripts/swarm_utils.py"


rule run_swarm:
    input:
        rules.format_swarm.output.fasta
    output:
        "results/swarm/{rundir}/swarm.txt"
    log:
        "logs/swarm/{rundir}/swarm.log"
    params:
        settings = config["swarm"]["settings"]
    threads: config["swarm"]["threads"]
    conda: "../envs/swarm.yml"
    shell:
        """
        swarm {params.settings} {input} -o {output} > {log} 2>&1
        """

rule swarm:
    input:
        expand("results/swarm/{rundir}/swarm.txt", rundir = config["rundir"])