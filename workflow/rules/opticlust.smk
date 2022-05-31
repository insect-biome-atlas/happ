rule mothur_sum:
    input:
        counts = expand("data/{rundir}/asv_counts.tsv", rundir = config["rundir"]),
        fasta = expand("data/{rundir}/asv_seqs.fasta", rundir=config["rundir"])
    output:
        counts = "results/opticlust/{rundir}/counts.tsv",
        fasta = "results/opticlust/{rundir}/asv_seqs.fasta"
    script:
        "../scripts/opticlust_utils.py"

rule zipfile:
    input:
        "{f}"
    output:
        "{f}.gz"
    shell:
        """
        gzip {input}
        """

rule mothur_align:
    input:
        fasta = "results/opticlust/{rundir}/asv_seqs.fasta.gz"
    output:
        dist = "results/opticlust/{rundir}/asv_seqs.dist.gz"
    log:
        log = "logs/opticlust/{rundir}/mothur_align.log",
        err = "logs/opticlust/{rundir}/mothur_align.err"
    conda:
        "../envs/opticlust.yml"
    params:
        indir = lambda wildcards, input: os.path.dirname(input.fasta[0]),
        tmpdir = "$TMPDIR/opticlust/{rundir}",
        fasta = "$TMPDIR/opticlust/{rundir}/asv_seqs.fasta"
    threads: config["opticlust"]["threads"]
    shell:
        """
        mkdir -p {params.tmpdir}
        gunzip -c {input.fasta} > {params.fasta}
        mothur "#set.dir(output={params.tmpdir});set.logfile(name={log.log}); pairwise.seqs(fasta={params.fasta}, processors={threads})" >{log.err} 2>&1
        gzip {params.tmpdir}/asv_seqs.dist
        mv {params.tmpdir}/asv_seqs.dist.gz {output.dist}
        rm -rf {params.tmpdir}
        """

rule vsearch_align:
    input:
        fasta = "results/opticlust/{rundir}/asv_seqs.fasta.gz"
    output:
        sim = "results/opticlust/{rundir}/asv_seqs.dist.gz"
    log:
        "logs/opticlust/{rundir}/vsearch_align.log"
    params:
        sim = "$TMPDIR/opticlust/{rundir}/asv_seqs.dist",
        fasta = "$TMPDIR/opticlust/{rundir}/asv_seqs.fasta",
        tmpdir = "$TMPDIR/opticlust/{rundir}"
    threads: config["opticlust"]["vsearch"]["threads"]
    conda:
        "../envs/opticlust.yml"
    shell:
        """
        mkdir -p {params.tmpdir}
        gunzip -c {input.fasta} > {params.fasta}
        vsearch --allpairs_global {params.fasta} --acceptall --threads {threads}\
            --userout {params.sim} --userfields "query+target+id" >{log} 2>&1
        gzip {params.sim}
        mv {params.sim}.gz {output.sim} 
        """

if config["opticlust"]["aligner"] == "vsearch":
    ruleorder: vsearch_align > mothur_align
    config["opticlust"]["sim"] = "sim=true,"
else:
    ruleorder: mothur_align > vsearch_align
    config["opticlust"]["sim"]= ""

rule run_opticlust:
    input:
        dist = "results/opticlust/{rundir}/asv_seqs.dist.gz",
        counts = "results/opticlust/{rundir}/counts.tsv"
    output:
        expand("results/opticlust/{{rundir}}/asv_seqs.opti_mcc.{suff}",
            suff = ["list", "sensspec", "steps"])
    log:
        log = "logs/opticlust/{rundir}/opticlust.log",
        err = "logs/opticlust/{rundir}/opticlust.err",
    params:
        dist = "$TMPDIR/opticlust/{rundir}/asv_seqs.dist",
        counts = "$TMPDIR/opticlust/{rundir}/counts.tsv",
        tmpdir = "$TMPDIR/opticlust/{rundir}",
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
        sim = config["opticlust"]["sim"]
    conda:
        "../envs/opticlust.yml"
    threads: config["opticlust"]["threads"]
    shell:
        """
        mkdir -p {params.tmpdir}
        gunzip -c {input.dist} > {params.dist} 
        cp {input.counts} {params.counts}
        mothur "#set.dir(output={params.tmpdir});set.logfile(name={log.log});cluster(column={params.dist}, count={params.counts}, {params.sim} method=opti)" >{log.err} 2>&1
        rm {params.dist} {params.counts}
        mv {params.tmpdir}/* {params.outdir}
        rm -rf {params.tmpdir}
        """

rule opticlust:
    input:
        expand("results/opticlust/{rundir}/asv_seqs.opti_mcc.{suff}",
            suff = ["list", "sensspec", "steps"], rundir = config["rundir"])
