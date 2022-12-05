localrules:
    format_swarm,
    swarm2tab,
    swarm,


rule format_swarm:
    input:
        fasta="results/common/{rundir}/{tax}/asv_seqs.fasta.gz",
        counts="results/common/{rundir}/{tax}/total_counts.tsv",
    output:
        fasta="results/swarm/{rundir}/{tax}/reformat.fasta.gz",
    params:
        tmpdir="$TMPDIR/{rundir}_{tax}_format_swarm",
        fasta="$TMPDIR/{rundir}_{tax}_format_swarm/reformat.fasta.gz",
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
        rules.format_swarm.output.fasta,
    output:
        expand(
            "results/swarm/{{rundir}}/{{tax}}/{{run_name}}/{f}",
            f=["swarm_table.tsv", "swarm.txt"],
        ),
    log:
        "logs/swarm/{rundir}/{tax}/swarm.{run_name}.log",
    params:
        fastidious=check_swarm_options(opt="fastidious"),
        differences=config["swarm"]["differences"],
        boundary=f"-b {config['swarm']['boundary']}"
        if config["swarm"]["boundary"] > 0
        else "",
        no_otu_breaking=check_swarm_options(opt="no-otu-breaking"),
        match_reward=f"-m {config['swarm']['match-reward']}"
        if config["swarm"]["differences"] > 1
        else "",
        mismatch_penalty=f"-p {config['swarm']['mismatch-penalty']}"
        if config["swarm"]["differences"] > 1
        else "",
        gap_opening_penalty=f"-g {config['swarm']['gap-opening-penalty']}"
        if config["swarm"]["differences"] > 1
        else "",
        gap_extension_penalty=f"-e {config['swarm']['gap-extension-penalty']}"
        if config["swarm"]["differences"] > 1
        else "",
        tmpdir="$TMPDIR/swarm/{rundir}/{tax}",
        fasta="$TMPDIR/swarm/{rundir}/{tax}/reformat.fasta",
        txt="$TMPDIR/swarm/{rundir}/{tax}/swarm.txt",
        tsv="$TMPDIR/swarm/{rundir}/{tax}/swarm_table.tsv",
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
    threads: config["swarm"]["threads"]
    conda:
        "../envs/swarm.yml"
    resources:
        runtime=60 * 24,
    shell:
        """
        exec &>{log}
        mkdir -p {params.tmpdir}
        gunzip -c {input} > {params.fasta}
        #swarm -d 0 -w {params.tmpdir}/derep.fasta -o /dev/null {params.fasta} 
        swarm {params.fastidious} {params.no_otu_breaking} -d {params.differences} {params.boundary} \
            {params.match_reward} {params.mismatch_penalty} {params.gap_opening_penalty} {params.gap_extension_penalty} \
            {params.fasta} -o {params.txt} -i {params.tsv} -t {threads}
        mv {params.txt} {params.outdir}
        mv {params.tsv} {params.outdir}
        rm -rf {params.tmpdir}
        """


rule swarm2tab:
    input:
        "results/swarm/{rundir}/{tax}/{run_name}/swarm.txt",
    output:
        "results/swarm/{rundir}/{tax}/{run_name}/asv_clusters.tsv",
    params:
        tmpdir="$TMPDIR/swarm/{rundir}/{tax}",
        out="$TMPDIR/swarm/{rundir}/{tax}/asv_clusters.tsv",
    script:
        "../scripts/swarm_utils.py"


rule swarm:
    input:
        expand(
            "results/swarm/{rundir}/{tax}/{run_name}/asv_clusters.tsv",
            rundir=config["rundir"],
            tax=taxa,
            run_name=config["run_name"],
        ),
