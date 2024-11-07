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
    conda: "../envs/vsearch.yml"
    container: "docker://quay.io/biocontainers/vsearch:2.29.1--h6a68c12_0"
    resources:
        runtime = 30,
    threads: 1
    shell:
        """
        vsearch --sintax {input.qry} --sintax_cutoff {params.cutoff} --randseed {params.seed} --db {input.db} --tabbedout {output} --threads 1 >{log} 2>&1
        """

def get_sintax_files(wildcards):
    checkpoint_output = checkpoints.split_input.get(**wildcards).output[0]
    return expand("results/taxonomy/sintax/{rundir}/splits/{split}.tab",
                    rundir=config["rundir"], 
                    split=glob_wildcards(os.path.join(checkpoint_output, "stdin.part_{split}.fasta")).split)

rule aggregate_sintax:
    """
    Concatenates the sintax output files into a single file
    """
    output:
        "results/taxonomy/sintax/{rundir}/sintax.tab"
    input:
        get_sintax_files,
    shell:
        """
        cat {input} > {output}
        """

rule parse_sintax:
    """
    Parses the sintax output file into a tsv file
    """
    output:
        "results/taxonomy/sintax/{rundir}/taxonomy.tsv"
    input:
        rules.aggregate_sintax.output
    log:
        "logs/sintax/{rundir}/parse_sintax.log"
    params:
        src=workflow.source_path("../scripts/sintax_tsv.py"),
        ranks=config["sintax"]["ranks"]
    shell:
        """
        python {params.src} -i {input} -o {output} -r {params.ranks} > {log} 2>&1
        """