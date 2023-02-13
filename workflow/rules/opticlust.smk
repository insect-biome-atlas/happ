localrules:
    opticlust,
    opticlust2tab,
    reformat_distmat,


rule mothur_align:
    input:
        fasta="results/common/{rundir}/{algo}/{tax}/asv_seqs.fasta.gz",
    output:
        dist="results/opticlust/{rundir}/{algo}/{tax}/asv_seqs.dist.gz",
    log:
        log="logs/opticlust/{rundir}/{algo}/{tax}/mothur_align.log",
        err="logs/opticlust/{rundir}/{algo}/{tax}/mothur_align.err",
    conda:
        "../envs/opticlust.yml"
    params:
        indir=lambda wildcards, input: os.path.dirname(input.fasta[0]),
        tmpdir="$TMPDIR/opticlust/{rundir}/{algo}/{tax}",
        fasta="$TMPDIR/opticlust/{rundir}/{algo}/{tax}/asv_seqs.fasta",
    threads: config["opticlust"]["threads"]
    resources:
        runtime=60 * 24 * 10,
    shell:
        """
        mkdir -p {params.tmpdir}
        gunzip -c {input.fasta} > {params.fasta}
        mothur "#set.dir(output={params.tmpdir});set.logfile(name={log.log}); pairwise.seqs(fasta={params.fasta}, processors={threads})" >{log.err} 2>&1
        gzip {params.tmpdir}/asv_seqs.dist
        mv {params.tmpdir}/asv_seqs.dist.gz {output.dist}
        rm -rf {params.tmpdir}
        """


def opticlust_input(wildcards):
    if config["opticlust"]["aligner"] == "vsearch":
        return f"results/vsearch/{wildcards.rundir}/{wildcards.algo}/{wildcards.tax}/asv_seqs.dist.reformat.gz"
    else:
        return f"results/opticlust/{wildcards.rundir}/{wildcards.algo}/{wildcards.tax}/asv_seqs.dist.gz"


rule reformat_distmat:
    input:
        "results/vsearch/{rundir}/{algo}/{tax}/asv_seqs.dist.gz",
    output:
        out="results/vsearch/{rundir}/{algo}/{tax}/asv_seqs.dist.reformat.gz",
    params:
        out="$TMPDIR/{rundir}_{algo}_{tax}_reformat_distmat/asv_seqs.dist.reformat.gz",
        tmpdir="$TMPDIR/{rundir}_{algo}_{tax}_reformat_distmat",
    script:
        "../scripts/opticlust_utils.py"


rule run_opticlust:
    """
    opticlust requires that all sequences in the counts file have abundance > 0
    """
    input:
        dist=opticlust_input,
        total_counts="results/common/{rundir}/{algo}/{tax}/total_counts.tsv",
    output:
        list="results/opticlust/{rundir}/{algo}/{tax}/{run_name}/asv_seqs.opti_mcc.list",
        sens="results/opticlust/{rundir}/{algo}/{tax}/{run_name}/asv_seqs.opti_mcc.sensspec",
        step="results/opticlust/{rundir}/{algo}/{tax}/{run_name}/asv_seqs.opti_mcc.steps",
    log:
        log="logs/opticlust/{rundir}/{algo}/{tax}/{run_name}/opticlust.log",
        err="logs/opticlust/{rundir}/{algo}/{tax}/{run_name}/opticlust.err",
    #shadow:
    #    "full"
    params:
        dist="$TMPDIR/opticlust.{rundir}.{algo}.{tax}/asv_seqs.dist",
        counts="$TMPDIR/opticlust.{rundir}.{algo}.{tax}/counts.tsv",
        tmpdir="$TMPDIR/opticlust.{rundir}.{algo}.{tax}",
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        delta=config["opticlust"]["delta"],
        cutoff=config["opticlust"]["cutoff"],
        initialize=config["opticlust"]["initialize"],
        precision=config["opticlust"]["precision"],
        src="workflow/scripts/run_opticlust.py",
    conda:
        "../envs/opticlust.yml"
    threads: config["opticlust"]["threads"]
    resources:
        runtime=60 * 24,
    shell:
        """
        mkdir -p {params.tmpdir} {params.outdir}
        gunzip -c {input.dist} > {params.dist} 
        cp {input.total_counts} {params.counts}
        python {params.src} {params.dist} {params.counts} --outdir {params.outdir} \
            --delta {params.delta} --initialize {params.initialize} --cutoff {params.cutoff} \
            --precision {params.precision} > {log.log}
        rm -rf {params.tmpdir}
        rm -f mothur.*.logfile
        """


rule opticlust2tab:
    """
    Generate a membership style table of clusters
    """
    input:
        rules.run_opticlust.output.list,
    output:
        "results/opticlust/{rundir}/{algo}/{tax}/{run_name}/asv_clusters.tsv",
    params:
        tmpdir="$TMPDIR/opticlust/{rundir}/{algo}/{tax}",
        out="$TMPDIR/opticlust/{rundir}/{algo}/{tax}/asv_clusters.tsv",
    script:
        "../scripts/opticlust_utils.py"


rule opticlust:
    input:
        expand(
            "results/opticlust/{rundir}/{algo}/{tax}/{run_name}/asv_clusters.tsv",
            rundir=config["rundir"],
            algo=config["chimera_algorithm"],
            tax=taxa,
            run_name=config["run_name"],
        ),
