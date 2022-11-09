localrules:
    opticlust,
    opticlust2tab,
    reformat_distmat,


rule mothur_align:
    input:
        fasta="results/common/{rundir}/{tax}/asv_seqs.fasta.gz",
    output:
        dist="results/opticlust/{rundir}/{tax}/asv_seqs.dist.gz",
    log:
        log="logs/opticlust/{rundir}/{tax}/mothur_align.log",
        err="logs/opticlust/{rundir}/{tax}/mothur_align.err",
    conda:
        "../envs/opticlust.yml"
    params:
        indir=lambda wildcards, input: os.path.dirname(input.fasta[0]),
        tmpdir="$TMPDIR/opticlust/{rundir}/{tax}",
        fasta="$TMPDIR/opticlust/{rundir}/{tax}/asv_seqs.fasta",
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
        return f"results/vsearch/{wildcards.rundir}/{wildcards.tax}/asv_seqs.dist.reformat.gz"
    else:
        return f"results/opticlust/{wildcards.rundir}/{wildcards.tax}/asv_seqs.dist.gz"


rule reformat_distmat:
    input:
        "results/vsearch/{rundir}/{tax}/asv_seqs.dist.gz",
    output:
        out="results/vsearch/{rundir}/{tax}/asv_seqs.dist.reformat.gz",
    params:
        out="$TMPDIR/{rundir}_{tax}_reformat_distmat/asv_seqs.dist.reformat.gz",
        tmpdir="$TMPDIR/{rundir}_{tax}_reformat_distmat",
    script:
        "../scripts/opticlust_utils.py"


rule run_opticlust:
    """
    opticlust requires that all sequences in the counts file have abundance > 0
    """
    input:
        dist=opticlust_input,
        total_counts="results/common/{rundir}/{tax}/total_counts.tsv",
    output:
        list="results/opticlust/{rundir}/{tax}/{run_name}/asv_seqs.opti_mcc.list",
        sens="results/opticlust/{rundir}/{tax}/{run_name}/asv_seqs.opti_mcc.sensspec",
        step="results/opticlust/{rundir}/{tax}/{run_name}/asv_seqs.opti_mcc.steps",
    log:
        log="logs/opticlust/{rundir}/{tax}/{run_name}/opticlust.log",
        err="logs/opticlust/{rundir}/{tax}/{run_name}/opticlust.err",
    shadow:
        "minimal"
    params:
        dist="$TMPDIR/opticlust/{rundir}/{tax}/asv_seqs.dist",
        counts="$TMPDIR/opticlust/{rundir}/{tax}/counts.tsv",
        tmpdir="$TMPDIR/opticlust/{rundir}/{tax}",
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        delta=config["opticlust"]["delta"],
        cutoff=config["opticlust"]["cutoff"],
        initialize=config["opticlust"]["initialize"],
        precision=config["opticlust"]["precision"],
    conda:
        "../envs/opticlust.yml"
    threads: config["opticlust"]["threads"]
    resources:
        runtime=60 * 24,
    shell:
        """
        # Check lines in input
        lines=$(gunzip -c {input.dist} | head | wc -l)
        if [ $lines -lt 1 ]; then
            echo "No results" > {output.list}
            cut -f1 {input.total_counts} | egrep -v "^Representative_Sequence" >> {output.list}
            touch {output.sens} {output.step} {log.log} {log.err}
        else
            mkdir -p {params.tmpdir}
            gunzip -c {input.dist} > {params.dist} 
            cp {input.total_counts} {params.counts}
            mothur "#set.dir(output={params.tmpdir});set.logfile(name={log.log});\
                cluster(column={params.dist}, count={params.counts}, \
                method=opti, delta={params.delta}, cutoff={params.cutoff}, initialize={params.initialize},\
                precision={params.precision})" >{log.err} 2>&1
            rm {params.dist} {params.counts}
            mv {params.tmpdir}/* {params.outdir}
            rm -rf {params.tmpdir}
        fi
        """


rule opticlust2tab:
    """
    Generate a membership style table of clusters
    """
    input:
        "results/opticlust/{rundir}/{tax}/{run_name}/asv_seqs.opti_mcc.list",
    output:
        "results/opticlust/{rundir}/{tax}/{run_name}/asv_clusters.tsv",
    params:
        tmpdir="$TMPDIR/opticlust/{rundir}/{tax}",
        out="$TMPDIR/opticlust/{rundir}/{tax}/asv_clusters.tsv",
    script:
        "../scripts/opticlust_utils.py"


rule opticlust:
    input:
        expand(
            "results/opticlust/{rundir}/{tax}/{run_name}/asv_clusters.tsv",
            rundir=config["rundir"],
            tax=taxa,
            run_name=config["opticlust"]["run_name"],
        ),
