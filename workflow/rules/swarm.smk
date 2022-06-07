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

def check_swarm_options(opt):
    if config["swarm"][opt]:
        return f"--{opt}"
    return ""

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
        fastidious = check_swarm_options(opt="fastidious"),
        differences = config["swarm"]["differences"],
        boundary = config["swarm"]["boundary"],
        no_otu_breaking = check_swarm_options(opt="no-otu-breaking"),
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
        swarm {params.fastidious} {params.no_otu_breaking} -d {params.differences} -b {params.boundary} \
            {params.fasta} -o {params.output} -t {threads} > {log} 2>&1
        mv {params.output} {output}
        rm -rf {params.tmpdir}
        """

rule swarm:
    input:
        expand("results/swarm/{rundir}/swarm.txt", rundir = config["rundir"])