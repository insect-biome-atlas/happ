import pandas as pd

def fetch_samples(f):
    r = pd.read_csv(f, sep="\t", nrows=1, index_col=0)
    return list(r.columns)

samples = fetch_samples(f=f"data/{config['rundir']}/asv_counts.tsv")

# Read fasta and counts file and split by sample

rule chimera:
    input:
        expand("data/{rundir}/per_sample/{sample}.sum.tsv",
            sample=samples,
            rundir=config["rundir"])

rule split_counts:
    output:
        expand("data/{{rundir}}/per_sample/{sample}.sum.tsv",
            sample=samples,)
    input:
        counts="data/{rundir}/asv_counts.tsv"
    log:
        "logs/chimeras/{rundir}/split_counts.log"
    params:
        src=srcdir("../scripts/split_counts.py"),
        outdir=lambda wildcards, output: os.path.dirname(output[0])
    shell:
        """
        python {params.src} {input.counts} {params.outdir} 2>{log}
        """

rule add_sums:
    output:
        fasta="data/{rundir}/per_sample/{sample}.fasta"
    input:
        sums="data/{rundir}/per_sample/{sample}.sum.tsv",
        fasta="data/{rundir}/asv_seqs.fasta"
    log:
        "logs/chimeras/{rundir}/{sample}.add-sums.log"
    params:
        src=srcdir("../scripts/add_size_to_fastaheader.py")
    shell:
        """
        python {params.src} {input.fasta} {input.sums} > {output.fasta} 2>{log}
        """

rule uchime_denovo:
    output:
        chim="results/chimera/{rundir}/{algo}/chimeras.fasta",
        nochim="results/chimera/{rundir}/{algo}/nonchimeras.fasta",
        border="results/chimera/{rundir}/{algo}/borderline.fasta",
        uchimeout="results/chimera/{rundir}/{algo}/uchimeout.txt",
    input:
        rules.add_sums.output.fasta
    log:
        "logs/chimeras/{rundir}.{algo}.log"
    conda:
        "../envs/vsearch.yml"
    threads: 4
        resources:
        runtime=60 * 24
    params:
        abskew=get_abskew,
        dn=config["vsearch"]["dn"],
        mindiffs=config["vsearch"]["mindiffs"],
        mindiv=config["vsearch"]["mindiv"],
        minh=config["vsearch"]["minh"],
    shell:
        """
        vsearch --threads {threads} --dn {params.dn} -- mindiffs {params.mindiffs} --mindiv {params.mindiv} --minh {params.minh} \
            {params.abskew} --chimeras {output.chim} --borderline {output.border} --nonchimeras {output.nochim} \
            --uchime_denovo {input.fasta} --uchimeout {output.uchimeout} >{log} 2>&1
        """