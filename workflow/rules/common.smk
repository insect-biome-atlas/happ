
localrules:
    filter_seqs

rule filter_seqs:
    input:
        counts = expand("data/{rundir}/asv_counts.tsv", rundir=config["rundir"]),
        fasta = expand("data/{rundir}/asv_seqs.fasta", rundir=config["rundir"])
    output:
        total_counts = "results/common/{rundir}/total_counts.tsv",
        counts = "results/common/{rundir}/asv_counts.tsv.gz",
        fasta = "results/common/{rundir}/asv_seqs.fasta.gz"
    script:
        "../scripts/common.py"

rule filter:
    input:
        expand("results/common/{rundir}/{f}",
            rundir = config["rundir"], f = ["total_counts.tsv", "asv_counts.tsv.gz", "asv_seqs.fasta.gz"])

rule vsearch_align:
    input:
        fasta = "results/common/{rundir}/asv_seqs.fasta.gz"
    output:
        dist = "results/vsearch/{rundir}/asv_seqs.dist.gz"
    log:
        "logs/vsearch/{rundir}/vsearch_align.log"
    params:
        dist = "$TMPDIR/vsearch/{rundir}/asv_seqs.dist",
        fasta="$TMPDIR/vsearch/{rundir}/asv_seqs.fasta",
        tmpdir="$TMPDIR/vsearch/{rundir}",
        id = config["vsearch"]["id"],
        iddef = config["vsearch"]["iddef"],
        query_cov = config["vsearch"]["query_cov"]
    threads: config["vsearch"]["threads"]
    conda:
        "../envs/opticlust.yml"
    resources:
        runtime = 60 * 24
    shell:
        """
        mkdir -p {params.tmpdir}
        gunzip -c {input.fasta} > {params.fasta}
        vsearch --usearch_global {params.fasta} --db {params.fasta} --self \
            --userout {params.dist} -userfields query+target+id --maxaccepts 0 --maxrejects 0 \
            --id {params.id} --iddef {params.iddef}  --query_cov {params.query_cov} --threads {threads} > {log} 2>&1
        gzip {params.dist}
        mv {params.dist}.gz {output.dist} 
        """