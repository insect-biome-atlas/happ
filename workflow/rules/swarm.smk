localrules:
    format_swarm,
    swarm2tab,
    swarm,


rule format_swarm:
    message: "Formatting files for SWARM for sequences in {wildcards.tax}"
    input:
        fasta="results/common/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/asv_seqs.fasta.gz",
        counts="results/common/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/total_counts.tsv",
    output:
        fasta="results/swarm/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/reformat.fasta.gz",
        derep=touch(
            "results/swarm/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/derep.txt"
        ),
    params:
        tmpdir=os.path.expandvars("$TMPDIR/{rundir}_{chimera_run}_{chimdir}_{rank}_{tax}_format_swarm"),
        fasta=os.path.expandvars("$TMPDIR/{rundir}_{chimera_run}_{chimdir}_{rank}_{tax}_format_swarm/reformat.fasta.gz"),
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
    message: "Running SWARM clustering on {wildcards.tax} sequences"
    input:
        rules.format_swarm.output.fasta,
    output:
        txt="results/swarm/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/{run_name}/swarm.txt",
        tsv="results/swarm/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/{run_name}/swarm_table.tsv",
    log:
        "logs/swarm/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/swarm.{run_name}.log",
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
        tmpdir="$TMPDIR/swarm/{rundir}/{chimera_run}/{chimdir}/{rank}/{tax}",
        fasta="$TMPDIR/swarm/{rundir}/{chimera_run}/{chimdir}/{rank}/{tax}/reformat.fasta",
        txt="$TMPDIR/swarm/{rundir}/{chimera_run}/{chimdir}/{rank}/{tax}/swarm.txt",
        tsv="$TMPDIR/swarm/{rundir}/{chimera_run}/{chimdir}/{rank}/{tax}/swarm_table.tsv",
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
    threads: 10
    conda: config["swarm-env"]
    container: "docker://quay.io/biocontainers/swarm:3.1.5--h4ac6f70_1"
    resources:
        runtime=60 * 24,
    shell:
        """
        exec &>{log}
        mkdir -p {params.tmpdir}
        gunzip -c {input} > {params.fasta}
        swarm {params.fastidious} {params.no_otu_breaking} -d {params.differences} {params.boundary} \
            {params.match_reward} {params.mismatch_penalty} {params.gap_opening_penalty} {params.gap_extension_penalty} \
            {params.fasta} -o {params.txt} -i {params.tsv} -t {threads}
        mv {params.txt} {params.outdir}
        mv {params.tsv} {params.outdir}
        rm -rf {params.tmpdir}
        """


rule swarm2tab:
    message: "Generating SWARM cluster file for {wildcards.tax}"
    input:
        rules.run_swarm.output.txt,
        rules.format_swarm.output.derep,
    output:
        "results/swarm/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/{run_name}/asv_clusters.tsv",
    params:
        tmpdir=os.path.expandvars("$TMPDIR/swarm/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}"),
        out=os.path.expandvars("$TMPDIR/swarm/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/asv_clusters.tsv"),
    script:
        "../scripts/swarm_utils.py"

def get_swarm_files(wildcards):
    checkpoint_dir = checkpoints.filter_seqs.get(
        rundir=config["rundir"], 
        chimera_run=config["chimera"]["run_name"], 
        chimdir=config["chimdir"], 
        rank=config["split_rank"]
        ).output[0]
    files = expand("results/swarm/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/{run_name}/asv_clusters.tsv",
        rundir=config["rundir"],
        chimera_run=config["chimera"]["run_name"],
        chimdir=config["chimdir"],
        rank=config["split_rank"],
        tax=glob_wildcards(os.path.join(checkpoint_dir, "{tax}", "asv_seqs.fasta.gz")).tax,
        run_name=config["run_name"]
        )
    return files

rule swarm:
    input:
        get_swarm_files,
        
