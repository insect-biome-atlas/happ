localrules:
    parse_sintax,
    aggregate_sintax,

rule sintax:
    """
    Runs sintax on one split of the query fasta file
    """
    output:
        temp("results/taxonomy/sintax/{rundir}/splits/{split}.tab")
    input:
        db=config["sintax"]["ref"],
        qry="results/common/{rundir}/splits/stdin.part_{split}.fasta"
    log:
        "logs/sintax/{rundir}/{split}.log"
    params:
        seed=config["sintax"]["randseed"],
        cutoff=config["sintax"]["cutoff"]
    conda: config["vsearch-env"]
    container: "docker://quay.io/biocontainers/vsearch:2.29.1--h6a68c12_0"
    resources:
        runtime = 30,
    threads: 1
    shell:
        """
        vsearch --sintax {input.qry} --sintax_cutoff {params.cutoff} --randseed {params.seed} --db {input.db} --tabbedout {output} --threads 1 >{log} 2>&1
        """

rule parse_sintax:
    """
    Parses the sintax output file into a tsv file
    """
    output:
        tsv=temp("results/taxonomy/sintax/{rundir}/splits/{split}.tsv"),
        conf=temp("results/taxonomy/sintax/{rundir}/splits/{split}_conf.tsv"),
    input:
        rules.sintax.output
    log:
        "logs/sintax/{rundir}/splits/parse_sintax.{split}.log"
    params:
        src=workflow.source_path("../scripts/sintax2tsv.py"),
        ranks=config["sintax"]["ranks"],
        cutoff=config["sintax"]["cutoff"]
    shell:
        """
        python {params.src} {input} {output.tsv} --conf_out {output.conf} -c {params.cutoff} -r {params.ranks} > {log} 2>&1
        """

def get_sintax_files(wildcards):
    checkpoint_output = checkpoints.split_input.get(**wildcards).output[0]
    d={}
    d["tsv"] = expand("results/taxonomy/sintax/{rundir}/splits/{split}.tsv",
                    rundir=config["rundir"], 
                    split=glob_wildcards(os.path.join(checkpoint_output, "stdin.part_{split}.fasta")).split)
    d["conf"] = [f.replace(".tsv", "_conf.tsv") for f in d["tsv"]]
    d["tab"] = [f.replace(".tsv", ".tab") for f in d["tsv"]]
    return d

rule aggregate_sintax:
    """
    Concatenates the sintax output files
    """
    output:
        tab="results/taxonomy/sintax/{rundir}/sintax.tab",
        tsv="results/taxonomy/sintax/{rundir}/taxonomy.tsv",
        conf="results/taxonomy/sintax/{rundir}/confidence.tsv"
    input:
        unpack(get_sintax_files)
    run:
        concat_files(input.tab).to_csv(output.tab, sep="\t")
        concat_files(input.tsv).to_csv(output.tsv, sep="\t")
        concat_files(input.conf).to_csv(output.conf, sep="\t")


