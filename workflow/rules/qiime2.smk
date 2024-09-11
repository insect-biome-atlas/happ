localrules:
    sintax2qiime_input,
    qiime2_import_qry_seqs,
    qiime2_import_ref_seqs,
    qiime2_import_taxonomy,
    qiime2_export,
    split_qiime_input,
    parse_qiime

wildcard_constraints:
    classifier = "vsearch|sklearn",

rule run_qiime2_vsearch:
    input:
        expand("results/taxonomy/vsearch/{rundir}/taxonomy.tsv",
            rundir = config["rundir"])

rule run_qiime2_sklearn:
    input:
        expand("results/taxonomy/sklearn/{rundir}/taxonomy.tsv",
            rundir = config["rundir"])

rule sintax2qiime_input:
    """
    Reformats a sintax database fasta to qiime2 input format
    Output should look like this:
    seqid1	k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Legionellales; f__Legionellaceae; g__Legionella; s__
    seqid2	k__Bacteria; p__Bacteroidetes; c__Flavobacteriia; o__Flavobacteriales; f__Flavobacteriaceae; g__Flavobacterium; s__

    fasta output should be:
    >seqid1
    ATGCGGGCTAGAGTAGCGAT...
    """
    input:
        fasta=config["sintax"]["ref"]
    output:
        tsv="resources/qiime2/sintax2qiime/taxonomy.tsv",
        fasta="resources/qiime2/sintax2qiime/seqs.fasta",
    log:
        "resources/qiime2/sintax2qiime/log"
    shell:
        """
        python workflow/scripts/sintax2qiime_input.py {input.fasta} {output.tsv} {output.fasta} > {log} 2>&1
        """

def qiime2_ref_seqs(wildcards):
    """
    Returns the path to the QIIME2 reference sequences if it exists,
    otherwise return reformatted sintax database
    """
    if os.path.exists(config["qiime2"]["ref"]):
        return config["qiime2"]["ref"]
    elif os.path.exists(config["sintax"]["ref"]):
        return "resources/qiime2/sintax2qiime/seqs.fasta"


rule qiime2_import_ref_seqs:
    """
    Imports the reference sequences into a QIIME2 artifact
    """
    output:
        "resources/qiime2/seqs.qza"
    input:
        qiime2_ref_seqs,
    log:
        "logs/qiime2_import_ref_seqs.log"
    container:
        "docker://quay.io/qiime2/core:2023.9" # may have to be built as part of a SLURM job on Uppmax
    threads: 1
    resources:
        runtime = 60
    shell:
        """
        qiime tools import --type 'FeatureData[Sequence]' --input-path {input} --output-path {output} > {log} 2>&1
        """
        
rule qiime2_import_qry_seqs:
    """
    Imports the query sequences into a QIIME2 artifact
    """
    output:
        temp("results/taxonomy/qiime2/{rundir}/splits/{split}/{split}.qza")
    input:
        "results/common/{rundir}/splits/stdin.part_{split}.fasta"
    log:
        "logs/qiime2_import_qry_seqs/{rundir}/{split}.log"
    container:
        "docker://quay.io/qiime2/core:2023.9" # may have to be built as part of a SLURM job on Uppmax
    threads: 1
    resources:
        runtime = 60
    shell:
        """
        qiime tools import --type 'FeatureData[Sequence]' --input-path {input} --output-path {output} > {log} 2>&1
        """

def qiime2_taxonomy(wildcards):
    """
    Returns the path to the QIIME2 taxonomy file if it exists,
    otherwise return reformatted sintax database
    """
    if os.path.exists(config["qiime2"]["taxfile"]):
        return config["qiime2"]["taxfile"]
    elif os.path.exists(config["sintax"]["ref"]):
        return "resources/qiime2/sintax2qiime/taxonomy.tsv"

rule qiime2_import_taxonomy:
    """
    Imports the taxonomy file into a QIIME2 artifact
    """
    output:
        "resources/qiime2/taxonomy.qza"
    input:
        qiime2_taxonomy,
    log:
        "logs/qiime2_import_taxonomy.log"
    container:
        "docker://quay.io/qiime2/core:2023.9" # may have to be built as part of a SLURM job on Uppmax
    threads: 1
    resources:
        runtime = 60
    shell:
        """
        qiime tools import --type 'FeatureData[Taxonomy]' --input-format TSVTaxonomyFormat \
            --input-path {input} --output-path {output} > {log} 2>&1
        """

rule qiime2_train:
    """
    Trains a naive bayes classifier on the reference sequences and taxonomy
    """
    output:
        "resources/qiime2/classifier.qza"
    input:
        tax=rules.qiime2_import_taxonomy.output[0],
        seq=rules.qiime2_import_ref_seqs.output[0]
    log:
        "logs/qiime2_train.log"
    container:
        "docker://quay.io/qiime2/core:2023.9" # may have to be built as part of a SLURM job on Uppmax
    resources:
        runtime=60 * 24 * 10,
    threads: 20
    shell:
        """
        qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads {input.seq} --i-reference-taxonomy {input.tax} \
            --o-classifier {output} > {log} 2>&1
        """

rule qiime2_classify_sklearn:
    """
    Classifies the query sequences using a naive bayes classifier
    """
    output:
        "results/taxonomy/sklearn/{rundir}/splits/{split}/taxonomy.qza"
    input:
        classifier=rules.qiime2_train.output[0],
        qry=rules.qiime2_import_qry_seqs.output[0]
    log:
        "logs/qiime2_classify_sklearn/{rundir}/{split}.log"
    threads: 20
    container:
        "docker://quay.io/qiime2/core:2023.9" # may have to be built as part of a SLURM job on Uppmax
    resources:
        runtime = 60 * 10,
    shell:
        """
        qiime feature-classifier classify-sklearn --i-classifier {input.classifier} --i-reads {input.qry} --o-classification {output} > {log} 2>&1
        """

rule qiime2_classify_vsearch:
    """
    Classifies the query sequences using vsearch
    """
    output:
        vsearch="results/taxonomy/vsearch/{rundir}/splits/{split}/taxonomy.qza",
        hits="results/taxonomy/vsearch/{rundir}/splits/{split}/hits.qza",
    input:
        ref=rules.qiime2_import_ref_seqs.output[0],
        ref_tax=rules.qiime2_import_taxonomy.output[0],
        qry=rules.qiime2_import_qry_seqs.output[0]
    log:
        "logs/qiime2_classify_vsearch/{rundir}/{split}.log"
    threads: 20
    container:
        "docker://quay.io/qiime2/core:2023.9" # may have to be built as part of a SLURM job on Uppmax
    resources:
        runtime = 60 * 10,
    shell:
        """
        qiime feature-classifier classify-consensus-vsearch --i-reference-reads {input.ref} --i-query {input.qry} \
            --i-reference-taxonomy {input.ref_tax} --o-classification {output.vsearch} --o-search-results {output.hits} \
            --p-threads {threads} --verbose > {log} 2>&1
        """

rule qiime2_export:
    output:
        "results/taxonomy/{classifier}/{rundir}/splits/{split}/taxonomy.tsv"
    input:
        "results/taxonomy/{classifier}/{rundir}/splits/{split}/taxonomy.qza"
    log:
        "logs/qiime2_export/{classifier}/{rundir}/{split}.log"
    container: 
        "docker://quay.io/qiime2/core:2023.9" # may have to be built as part of a SLURM job on Uppmax
    threads: 1
    shell:
        """
        qiime tools export --input-path {input} --output-path {output[0]} --output-format TSVTaxonomyFormat > {log} 2>&1
        """

def get_classifier_files(wildcards):
    checkpoint_output = checkpoints.split_input.get(**wildcards).output[0]
    return expand("results/taxonomy/{classifier}/{rundir}/splits/{split}/taxonomy.tsv",
                    classifier=wildcards.classifier, rundir=wildcards.rundir,
                    split=glob_wildcards(os.path.join(checkpoint_output, "stdin.part_{split}.fasta")).split)

rule aggregate_qiime:
    """
    Concatenates the qiime output files into a single file
    """
    output:
        "results/taxonomy/{classifier}/{rundir}/taxonomy.raw.tsv"
    input:
        get_classifier_files,
    run:
        with open(output[0], "w") as out:
            for i, f in enumerate(input):
                with open(f, 'r') as infile:
                    for j, line in enumerate(infile):
                        if j == 0 and i == 0:
                            out.write(line)
                        elif j > 0:
                            out.write(line)

def add_unassigned(row):
    _row = row[0:-1]
    if sum(_row!="") == 0:
        return row.replace("", "unclassified")
    last_known = _row[_row!=''].index[-1]
    if last_known == _row.index[-1]:
        return row
    if _row[last_known].startswith("unclassified"):
        uncl_str = _row[last_known]
    else:
        uncl_str = f"unclassified.{_row[_row.index[_row!=''][-1]]}"
    row[row==""] =  uncl_str
    return row

rule parse_qiime:
    output:
        "results/taxonomy/{classifier}/{rundir}/taxonomy.tsv"
    input:
        rules.aggregate_qiime.output[0]
    run:
        import pandas as pd
        df = pd.read_csv(input[0], sep="\t", index_col=0)
        parsed = parse_qiime2(df)
        parsed = parsed.apply(add_unassigned, axis=1)
        parsed.to_csv(output[0], sep="\t")
    