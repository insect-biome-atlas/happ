localrules:
    parse_sintax,
    collate_sintax,
    extract_ASVs,
    split_sintax_input

checkpoint split_sintax_input:
    """
    Splits the sintax fasta file into 1000 chunks
    """
    output:
        directory("results/sintax/{rundir}/splits")
    input:
        qry="data/{rundir}/asv_seqs.fasta"
    log:
        "logs/split_sintax_input/{rundir}.split.log"
    params:
        outdir=lambda wildcards, output: output[0],
        size=500,
    resources:
        runtime = 60,
    threads: 1
    shell:
        """
        cat {input.qry} | seqkit split2 -O {params.outdir} -j {threads} -s {params.size} >{log} 2>&1
        """

rule sintax:
    """
    Runs sintax on one split of the query fasta file
    """
    output:
        temp("results/sintax/{rundir}/splits/{split}.tab")
    input:
        db=config["sintax"]["ref"]["fasta"],
        qry="results/sintax/{rundir}/splits/stdin.part_{split}.fasta"
    log:
        "logs/sintax/{rundir}/{split}.log"
    params:
        seed=config["sintax"]["randseed"],
        cutoff=config["sintax"]["cutoff"]
    conda: "../envs/sintax.yml"
    resources:
        runtime = 30,
    threads: 1
    shell:
        """
        vsearch --sintax {input.qry} --sintax_cutoff {params.cutoff} --randseed {params.seed} --db {input.db} --tabbedout {output} --threads 1 >{log} 2>&1
        """

def aggregate_sintax(wildcards):
    checkpoint_output = checkpoints.split_sintax_input.get(**wildcards).output[0]
    return expand("results/sintax/{rundir}/splits/{split}.tab",
                    rundir=config["rundir"], 
                    split=glob_wildcards(os.path.join(checkpoint_output, "stdin.part_{split}.fasta")).split)

rule collate_sintax:
    """
    Concatenates the sintax output files into a single file
    """
    output:
        "results/sintax/{rundir}/sintax.tab"
    input:
        aggregate_sintax,
    shell:
        """
        cat {input} > {output}
        """

rule parse_sintax:
    """
    Parses the sintax output file into a tsv file
    """
    output:
        "results/sintax/{rundir}/sintax.tsv"
    input:
        rules.collate_sintax.output
    log:
        "logs/sintax/{rundir}/parse_sintax.log"
    params:
        src=workflow.source_path("../scripts/sintax_tsv.py"),
        ranks=lambda wildcards: config["sintax"]["ranks"]
    shell:
        """
        python {params.src} -i {input} -o {output} -r {params.ranks} > {log} 2>&1
        """

rule extract_ASVs:
    """
    Extract ASVs matching reassign_rank == reassign_taxa but unclassified at child rank.
    For example, ASVs classified as reassign_taxa == 'Insecta' at reassign_rank == 'Class' but unclassified at 'Order'.
    """
    output:
        tsv="results/sintax/{rundir}/reassign/{rank}/unclassified.{taxa}.tsv",
        fasta="results/sintax/{rundir}/reassign/{rank}/unclassified.{taxa}.fasta"
    input:
        tsv=rules.parse_sintax.output[0],
        qry="data/{rundir}/asv_seqs.fasta"
    params:
        reassign_rank = lambda wildcards: wildcards.rank,
        reassign_taxa = lambda wildcards: wildcards.taxa,
    run:
        reassign_taxa = params.reassign_taxa
        reassign_rank = params.reassign_rank
        from Bio.SeqIO import parse, write as write_fasta
        import pandas as pd
        df = pd.read_csv(input.tsv, sep="\t", index_col=0)
        cols = df.columns.tolist()
        child_rank = cols[cols.index(reassign_rank)+1]
        taxdf = df.loc[df[child_rank]==f"unclassified.{reassign_taxa}"]
        seqs = []
        for record in parse(input.qry, "fasta"):
            if record.id in taxdf.index:
                seqs.append(record)
        with open(output.fasta, "w") as f:
            write_fasta(seqs, f, "fasta")
        taxdf.to_csv(output.tsv, sep="\t")