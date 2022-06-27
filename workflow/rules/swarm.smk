localrules:
    format_swarm,
    swarm

rule format_swarm:
    input:
        fasta = "results/common/{rundir}/{tax}/asv_seqs.fasta.gz",
        counts = "results/common/{rundir}/{tax}/total_counts.tsv"
    output:
        fasta = "results/swarm/{rundir}/{tax}/reformat.fasta.gz"
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
        expand("results/swarm/{{rundir}}/{{tax}}/{f}", f = ["swarm_table.tsv", "swarm.txt"])
    log:
        "logs/swarm/{rundir}/{tax}/swarm.log"
    params:
        fastidious = check_swarm_options(opt="fastidious"),
        differences = config["swarm"]["differences"],
        boundary = config["swarm"]["boundary"],
        no_otu_breaking = check_swarm_options(opt="no-otu-breaking"),
        tmpdir = "$TMPDIR/swarm/{rundir}/{tax}",
        fasta = "$TMPDIR/swarm/{rundir}/{tax}/reformat.fasta",
        txt = "$TMPDIR/swarm/{rundir}/{tax}/swarm.txt",
        tsv = "$TMPDIR/swarm/{rundir}/{tax}/swarm_table.tsv",
        outdir = lambda wildcards, output: os.path.dirname(output[0])
    threads: config["swarm"]["threads"]
    conda: "../envs/swarm.yml"
    resources:
        runtime = 60 * 24 * 10
    shell:
        """
        mkdir -p {params.tmpdir}
        gunzip -c {input} > {params.fasta}
        swarm {params.fastidious} {params.no_otu_breaking} -d {params.differences} -b {params.boundary} \
            {params.fasta} -o {params.txt} -i {params.tsv} -t {threads} > {log} 2>&1
        mv {params.txt} {params.outdir}
        mv {params.tsv} {params.outdir}
        rm -rf {params.tmpdir}
        """

rule swarm2tab:
    input:
        "results/swarm/{rundir}/{tax}/swarm_table.tsv"
    output:
        "results/swarm/{rundir}/{tax}/asv_clusters.tsv"
    params:
        tmpdir = "$TMPDIR/swarm/{rundir}/{tax}",
        out = "$TMPDIR/swarm/{rundir}/{tax}/asv_clusters.tsv"
    script:
        "../scripts/swarm_utils.py"

rule swarm:
    input:
        expand("results/swarm/{rundir}/{tax}/asv_clusters.tsv",
            rundir = config["rundir"], tax=taxa)