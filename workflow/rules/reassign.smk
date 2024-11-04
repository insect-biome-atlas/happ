localrules:
    update_taxonomy

rule update_taxonomy:
    output:
        "results/taxonomy/sintax_epang/{rundir}/{heur}/taxonomy.tsv"
    input:
        sintax="results/taxonomy/sintax/{rundir}/taxonomy.tsv",
        epang="results/taxonomy/epa-ng/{rundir}/assignments/{heur}/taxonomy.tsv"
    params:
        agree_rank = lambda wildcards: config["reassign"]["placeholder_rank"],
        agree_taxa = lambda wildcards: config["reassign"]["placeholder_taxa"],
        update_ranks = lambda wildcards: config["reassign"]["reassign_ranks"],
        downstream_ranks = lambda wildcards: config["reassign"]["downstream_ranks"],
        src=workflow.source_path("../scripts/update_taxonomy.py")
    log:
        "logs/update_taxonomy/{rundir}/{heur}.log"
    shell:
        """
        python {params.src} -b {input.sintax} -u {input.epang} \
            -a {params.agree_rank} -t {params.agree_taxa} -U {params.update_ranks} -d {params.downstream_ranks} -o {output} 2>{log}
        """