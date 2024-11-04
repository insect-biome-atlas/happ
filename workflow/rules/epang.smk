localrules:
    nexus2newick,
    extract_ref_taxonomy,
    nexus2fasta,
    split_aln,
    gappa2taxdf,
    collate_gappa_taxdf,
    filter_query_aln,
    add_missing

wildcard_constraints:
    heur="baseball-heur|dyn-heur|no-heur",

## target rules

""" def epa_ng_input(wildcards):
    input = []
    for ref in config["phylogeny"]["ref"].keys():
        if "query" in config["phylogeny"].keys():
            for query in config["phylogeny"]["query"].keys():
                for heur in config["phylogeny"]["epa-ng"]["heuristics"]:
                    input.append(f"results/epa-ng/{ref}/queries/{query}/{heur}/taxonomy.tsv")
        if "reassign_sintax" in config["phylogeny"]["ref"][ref].keys():
            rank = config["phylogeny"]["ref"][ref]["reassign_sintax"]["rank"]
            taxa = config["phylogeny"]["ref"][ref]["reassign_sintax"]["taxa"]
            for sintax_ref in config["sintax"]["ref"].keys():
                for sintax_query in config["sintax"]["query"].keys():
                    for tax in taxa:
                        new_query = f"{sintax_query}.{tax}.reassigned"
                        for heur in config["phylogeny"]["epa-ng"]["heuristics"]:
                            input.append(f"results/epa-ng/{ref}/queries/{new_query}/{heur}/taxonomy.tsv")
                            input.append(f"results/reassign/{sintax_ref}/queries/{sintax_query}/{ref}/{heur}/taxonomy.tsv")
    return input """

#rule run_epa_ng:
#    input:
#        epa_ng_input,

rule nexus2newick:
    """
    Converts a nexus tree to newick format
    """
    output:
        "resources/epa-ng/tree.nwk",
    input:
        config["epa-ng"]["tree"],
    log:
        "logs/nexus2newick.log",
    params:
        src=workflow.source_path("../scripts/nexus2newick.py"),
    shell:
        """
        python {params.src} {input} {output} >{log} 2>&1
        """

def ref_tree(wildcards):
    if config["epa-ng"]["tree_format"] == "nexus":
        return rules.nexus2newick.output[0]
    elif config["epa-ng"]["tree_format"] == "newick":
        return config["epa-ng"]["tree"]

def ref_msa(wildcards):
    if config["epa-ng"]["msa_format"] == "nexus":
        return rules.nexus2fasta.output[0]
    elif config["epa-ng"]["msa_format"] == "fasta":
        return config["epa-ng"]

rule extract_ref_taxonomy:
    output:
        "resources/epa-ng/taxon_file.tsv",
    input:
        ref_tree,
    log:
        "logs/epa-ng/extract_ref_taxonomy.log",
    params:
        config["epa-ng"]["tree_ranks"],
        src=workflow.source_path("../scripts/extract_ref_taxonomy.py"),
    shell:
        """
        python {params.src} {input} {output} --ranks {params.ranks} >{log} 2>&1
        """

rule nexus2fasta:
    output:
        "resources/epa-ng/ref_msa.fasta",
    input:
        config["epa-ng"]["msa"],
    log:
        "logs/epa-ng/nexus2fasta.log",
    params:
        src=workflow.source_path("../scripts/convertalign.py"),
    shell:
        """
        python {params.src} {input} nexus {output} fasta >{log} 2>&1
        """

def ref_msa(wildcards):
    if config["epa-ng"]["msa_format"] == "nexus":
        return rules.nexus2fasta.output[0]
    elif config["epa-ng"]["msa_format"] == "fasta":
        return config["epa-ng"]["msa"]

rule hmm_build:
    output:
        "resources/epa-ng/ref_msa.hmm",
    input:
        ref_msa,
    log:
        "logs/epa-ng/hmmbuild.log",
    conda:
        "../envs/epa-ng.yml"
    resources:
        runtime=60,
    shell:
        """
        hmmbuild {output} {input} > {log} 2>&1
        """

rule hmm_align:
    output:
        temp("results/taxonomy/epa-ng/{rundir}/hmmalign/splits/{split}/{split}.fasta"),
    input:
        hmm=rules.hmm_build.output,
        qry="results/common/{rundir}/splits/stdin.part_{split}.fasta",
        ref_msa=ref_msa
    log:
        "logs/epa-ng/{rundir}/hmmalign/{split}.log",
    conda:
        "../envs/epa-ng.yml"
    resources:
        runtime=5,
    params:
        tmpdir=lambda wildcards: f"$TMPDIR/{wildcards.rundir}.{wildcards.split}.hmm_align",
    shell:
        """
        mkdir -p {params.tmpdir}
        hmmalign --trim --mapali {input.ref_msa} --outformat afa -o {params.tmpdir}/output.msa {input.hmm} {input.qry} > {log} 2>&1
        cut -f1 -d ' ' {params.tmpdir}/output.msa > {output[0]}
        rm -r {params.tmpdir}
        """

rule split_aln:
    output:
        ref_msa=temp("results/taxonomy/epa-ng/{rundir}/hmmalign/splits/{split}/reference.fasta"),
        qry_msa=temp("results/taxonomy/epa-ng/{rundir}/hmmalign/splits/{split}/query.fasta"),
    input:
        ref_msa=ref_msa,
        msa=rules.hmm_align.output,
    log:
        "logs/epa-ng/{rundir}/split_aln/{split}.log",
    params:
        outdir=lambda wildcards, output: os.path.dirname(output.ref_msa),
    conda: "../envs/epa-ng.yml"
    shell:
        """
        epa-ng --redo --split {input.ref_msa} {input.msa} --outdir {params.outdir} > {log} 2>&1
        """

rule filter_query_aln:
    """
    Removes sequences in the alignment that have too few positions aligned
    """
    output:
        keep=temp("results/taxonomy/epa-ng/{rundir}/hmmalign/splits/{split}/query.filtered.fasta"),
        r=temp("results/taxonomy/epa-ng/{rundir}/hmmalign/splits/{split}/query.removed.txt")
    input:
        rules.split_aln.output.qry_msa,
    params:
        min_len=50 #TODO: make this a parameter
    run:
        from Bio.AlignIO import read
        from Bio.SeqIO import write
        alignment = read(input[0], "fasta")
        keep = []
        remove = []
        for record in alignment:
            seq =record.seq
            if len(seq.replace("-", "").replace(".", "")) >= params.min_len:
                keep.append(record)
            else:
                remove.append(record.id)
        write(keep, output.keep, "fasta")
        with open(output.r, "w") as f:
            for r in remove:
                f.write(f"{r}\n")                     

rule raxml_evaluate:
    output:
        temp("results/taxonomy/epa-ng/{rundir}/raxml-ng/splits/{split}/info.raxml.bestModel"),
    input:
        tree=ref_tree,
        msa=rules.split_aln.output.ref_msa,
    log:
        "logs/epa-ng/{rundir}/raxml_evaluate/{split}.log",
    conda:
        "../envs/epa-ng.yml"
    params:
        model=lambda wildcards: config["epa-ng"]["model"],
        prefix=lambda wildcards, output: os.path.dirname(output[0]) + "/info",
    threads: 2
    shell:
        """
        raxml-ng --extra thread-pin --redo --threads {threads} --evaluate --msa {input.msa} --tree {input.tree} --prefix {params.prefix} --model {params.model} >{log} 2>&1
        """

## epa-ng

def get_heuristic(wildcards):
    if wildcards.heur == "no-heur":
        return "--no-heur"
    elif wildcards.heur == "dyn-heur":
        return "--dyn-heur 0.99999"
    elif wildcards.heur == "baseball-heur":
        return "--baseball-heur"

rule epa_ng:
    output:
        temp("results/taxonomy/epa-ng/{rundir}/placements/splits/{split}/{heur}/epa-ng_result.jplace"),
    input:
        qry=rules.filter_query_aln.output[0],
        ref_msa=rules.split_aln.output.ref_msa,
        ref_tree=ref_tree,
        info=rules.raxml_evaluate.output[0],
    log:
        "logs/epa-ng/{rundir}/placements/{split}/{heur}/epa-ng.log",
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        heur=get_heuristic,
    conda: 
        "../envs/epa-ng.yml"
    threads: 20
    resources:
        runtime=60*24,
    shell:
        """
        epa-ng --redo -T {threads} --tree {input.ref_tree} --ref-msa {input.ref_msa} \
            --query {input.qry} --out-dir {params.outdir} {params.heur} --model {input.info} >{log} 2>&1
        mv {params.outdir}/epa_result.jplace {output[0]}
        """

## gappa

def ref_taxonomy(wildcards):
    if config["epa-ng"]["ref_taxonomy"]:
        return config["epa-ng"]["ref_taxonomy"]
    else:
        return rules.extract_ref_taxonomy.output


def get_dist_ratio(config):
    dist_ratio = config["epa-ng"]["gappa"]["distribution_ratio"]
    if dist_ratio == -1:
        return ""
    else:
        return f"--distribution-ratio {dist_ratio}"

rule gappa_assign:
    """
    Run gappa taxonomic assignment on placement file
    """
    output:
        temp("results/taxonomy/epa-ng/{rundir}/assignments/splits/{split}/{heur}/per_query.tsv"),
        temp("results/taxonomy/epa-ng/{rundir}/assignments/splits/{split}/{heur}/profile.tsv"),
        temp("results/taxonomy/epa-ng/{rundir}/assignments/splits/{split}/{heur}/labelled_tree.newick"),
    input:
        jplace=rules.epa_ng.output,
        taxonfile=ref_taxonomy,
    log:
        "logs/epa-ng/{rundir}/assignments/splits/{split}/{heur}/gappa_assign.log",
    params:
        ranks_string=lambda wildcards: "|".join(config["epa-ng"]["tree_ranks"]),
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        consensus_thresh=config["epa-ng"]["gappa"]["consensus_thresh"],
        distribution_ratio=get_dist_ratio(config),
    conda:
        "../envs/epa-ng.yml"
    threads: 4
    #resources:
        #runtime=20,
        #cpus_per_task=1,
        #tasks=4
    shell:
        """
        gappa examine assign --threads {threads} --out-dir {params.outdir} \
            --jplace-path {input.jplace} --taxon-file {input.taxonfile} \
            --ranks-string '{params.ranks_string}' --per-query-results \
            --consensus-thresh {params.consensus_thresh} {params.distribution_ratio} \
            --best-hit --allow-file-overwriting > {log} 2>&1
        """


rule gappa2taxdf:
    """
    Convert gappa output to a taxonomic dataframe
    """
    output:
        temp("results/taxonomy/epa-ng/{rundir}/assignments/splits/{split}/{heur}/taxonomy.tsv"),
    input:
        rules.gappa_assign.output[0],
    log:
        "logs/epa-ng/{rundir}/assignments/splits/{split}/{heur}/gappa2taxdf.log"
    params:
        src=workflow.source_path("../scripts/gappa2taxdf.py"),
        ranks=lambda wildcards: config["epa-ng"]["tree_ranks"],
    shell:
        """
        python {params.src} {input} {output} --ranks {params.ranks} >{log} 2>&1
        """

def aggregate_gappa(wildcards):
    checkpoint_output = checkpoints.split_input.get(**wildcards).output[0]
    return expand("results/taxonomy/epa-ng/{rundir}/assignments/splits/{split}/{heur}/taxonomy.tsv",
                    rundir=config["rundir"],
                    heur=wildcards.heur,
                    split=glob_wildcards(os.path.join(checkpoint_output, "stdin.part_{split}.fasta")).split)

rule collate_gappa_taxdf:
    output:
        "results/taxonomy/epa-ng/{rundir}/assignments/{heur}/taxonomy.raw.tsv"
    input:
        aggregate_gappa
    run:
        with open(output[0], "w") as out:
            for i, f in enumerate(input):
                with open(f, 'r') as infile:
                    for j, line in enumerate(infile):
                        if j == 0 and i == 0:
                            out.write(line)
                        elif j > 0:
                            out.write(line)

def aggregate_removed(wildcards):
    checkpoint_output = checkpoints.split_input.get(**wildcards).output[0]
    return expand("results/taxonomy/epa-ng/{rundir}/hmmalign/splits/{split}/query.removed.txt",
                    rundir=config["rundir"], 
                    split=glob_wildcards(os.path.join(checkpoint_output, "stdin.part_{split}.fasta")).split)

rule add_missing:
    """
    Adds ASVs to the taxonomy table that were filtered out during the alignment step
    """
    output:
        "results/taxonomy/epa-ng/{rundir}/assignments/{heur}/taxonomy.tsv"
    input:
        df=rules.collate_gappa_taxdf.output[0],
        asv=aggregate_removed,
    run:
        import pandas as pd
        df = pd.read_csv(input.df, sep="\t", index_col=0)
        for f in input.asv:
            with open(f, "r") as fhin:
                for line in fhin:
                    asv = line.rstrip()
                    if asv not in df.index:
                        df.loc[asv] = ["unclassified"] * len(df.columns)
        df.to_csv(output[0], sep="\t")