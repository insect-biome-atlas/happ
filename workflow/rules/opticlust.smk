localrules:
    opticlust,
    opticlust2tab,
    reformat_distmat,


rule mothur_align:
    input:
        fasta="results/common/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/asv_seqs.fasta.gz",
    output:
        dist="results/mothur/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/asv_seqs.dist.gz",
    log:
        log="logs/mothur/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/mothur_align.log",
        err="logs/mothur/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/mothur_align.err",
    conda: config["opticlust-env"]
    container: "docker://biocontainers/mothur:v1.41.21-1-deb_cv1"
    params:
        indir=lambda wildcards, input: os.path.dirname(input.fasta[0]),
        tmpdir="$TMPDIR/opticlust/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}",
        fasta="$TMPDIR/opticlust/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/asv_seqs.fasta",
    threads: 10
    resources:
        tasks = 10
    shell:
        """
        mkdir -p {params.tmpdir}
        gunzip -c {input.fasta} > {params.fasta}
        mothur "#set.dir(output={params.tmpdir});set.logfile(name={log.log}); pairwise.seqs(fasta={params.fasta}, processors={resources.tasks})" >{log.err} 2>&1
        gzip {params.tmpdir}/asv_seqs.dist
        mv {params.tmpdir}/asv_seqs.dist.gz {output.dist}
        rm -rf {params.tmpdir}
        """


def opticlust_input(wildcards):
    if config["opticlust"]["aligner"] == "vsearch":
        return f"results/vsearch/{wildcards.rundir}/{wildcards.chimera_run}/{wildcards.chimdir}/{wildcards.rank}/taxa/{wildcards.tax}/asv_seqs.dist.reformat.gz"
    else:
        return f"results/mothur/{wildcards.rundir}/{wildcards.chimera_run}/{wildcards.chimdir}/{wildcards.rank}/taxa/{wildcards.tax}/asv_seqs.dist.gz"


rule reformat_distmat:
    input:
        rules.vsearch_align.output.dist,
    output:
        out="results/vsearch/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/asv_seqs.dist.reformat.gz",
    params:
        out=os.path.expandvars("$TMPDIR/{rundir}_{chimera_run}_{chimdir}_{tax}_reformat_distmat/asv_seqs.dist.reformat.gz"),
        tmpdir=os.path.expandvars("$TMPDIR/{rundir}_{chimera_run}_{chimdir}_{tax}_reformat_distmat"),
    script:
        "../scripts/opticlust_utils.py"


rule run_opticlust:
    """
    opticlust requires that all sequences in the counts file have abundance > 0
    """
    input:
        dist=opticlust_input,
        counts="results/common/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/total_counts.tsv",
    output:
        list=touch("results/opticlust/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/{run_name}/asv_seqs.opti_mcc.list"),
        sens=touch("results/opticlust/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/{run_name}/asv_seqs.opti_mcc.sensspec"),
        step=touch("results/opticlust/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/{run_name}/asv_seqs.opti_mcc.steps"),
    log:
        "logs/opticlust/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/{run_name}/opticlust.log",
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        delta=config["opticlust"]["delta"],
        cutoff=config["opticlust"]["cutoff"],
        initialize=config["opticlust"]["initialize"],
        precision=config["opticlust"]["precision"],
        #src=workflow.source_path("../scripts/run_opticlust.py"),
    conda: config["opticlust-env"]
    container: "docker://biocontainers/mothur:v1.41.21-1-deb_cv1"
    shadow: "minimal"
    shell:
        """
        gunzip -c {input.dist} > asv_seqs.dist
        if [ -s asv_seqs.dist ]; then
            mothur "#set.dir(output={params.outdir});cluster(column=asv_seqs.dist, count={input.counts},method=opti, delta={params.delta}, cutoff={params.cutoff}, initialize={params.initialize},precision={params.precision})" > {log} 2>&1
        else
            echo "No results" > {output.list}
            cat {input.counts} | egrep -v "^Representative_Sequence" | cut -f1 >> {output.list}
        fi
        """


rule opticlust2tab:
    """
    Generate a membership style table of clusters
    """
    input:
        rules.run_opticlust.output.list,
    output:
        "results/opticlust/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/{run_name}/asv_clusters.tsv",
    params:
        tmpdir=os.path.expandvars("$TMPDIR/opticlust/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}"),
        out=os.path.expandvars("$TMPDIR/opticlust/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/asv_clusters.tsv"),
    script:
        "../scripts/opticlust_utils.py"

def get_opticlust_files(wildcards):
    checkpoint_dir = checkpoints.filter_seqs.get(
        rundir=config["rundir"], 
        chimera_run=config["chimera"]["run_name"], 
        chimdir=config["chimdir"], 
        rank=config["split_rank"]
        ).output[0]
    files = expand("results/opticlust/{rundir}/{chimera_run}/{chimdir}/{rank}/taxa/{tax}/{run_name}/asv_clusters.tsv",
        rundir=config["rundir"],
        chimera_run=config["chimera"]["run_name"],
        chimdir=config["chimdir"],
        rank=config["split_rank"],
        tax=glob_wildcards(os.path.join(checkpoint_dir, "{tax}", "asv_seqs.fasta.gz")).tax,
        run_name=config["run_name"]
        )
    return files

rule opticlust:
    input:
        get_opticlust_files,
