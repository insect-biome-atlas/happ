localrules:
    generate_counts_files,
    generate_taxa_seqs,
    generate_aa_seqs,
    generate_datasets,
    generate_evodistlists

rule generate_counts_files:
    """
    Generate the counts files for use with neeat filtering
    """
    output:
        cluster_counts="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/counts/cluster_counts.tsv",
        calibrated_counts="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/counts/calibrated_cluster_counts.tsv",
        tot_proportional_counts="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/counts/tot_proportional_cluster_counts.tsv.tsv",
        sample_proportional_counts="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/counts/sample_proportional_cluster_counts.tsv",
    input:
        counts="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_counts.tsv",
        taxonomy="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_taxonomy.tsv",
    log:
        "logs/generate_counts_files/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}.log"
    params:
        meta=config["metadata"]["metadata_file"] if os.path.exists(config["metadata"]["metadata_file"]) else False,
        sample_id_col=config["metadata"]["sample_id_col"],
        sample_type_col=config["metadata"]["sample_type_col"],
        sample_val=config["metadata"]["sample_val"],
        functions=workflow.source_path("../scripts/neeat/spikes_controls_fxns.R")
    conda:
        "../envs/r-env.yml"
    script: 
        "../scripts/neeat/generate_counts_files.R"

checkpoint generate_taxa_seqs:
    """
    Outputs a fasta file for each taxon (e.g. order) in the dataset. Skips taxa with 
    fewer than 2 ASVs.
    """
    output:
        directory("results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/fasta"),
    input:
        taxonomy="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_taxonomy.tsv",
        fasta="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_reps.fasta",
    shadow: "minimal"
    run:
        import subprocess
        import pandas as pd
        outdir=output[0]
        os.makedirs(outdir, exist_ok=True)
        # read in the taxonomy file
        # translate rank and columns to lowercase
        rank = (wildcards.noise_rank).lower()
        taxdf = pd.read_csv(input.taxonomy, sep="\t", index_col=0)
        taxdf.rename(columns = lambda x: x.lower(), inplace=True)
        # remove unclassified and ambiguous taxa
        taxdf = taxdf.loc[(~taxdf[rank].str.contains("_X+$"))&(~taxdf[rank].str.startswith("unclassified"))]
        # extract representative ASVs
        taxdf = taxdf.loc[taxdf["representative"]==1]
        # taxa is the unique set of values for split_rank
        taxa = taxdf[rank].unique()
        # iterate over taxa
        for tax in taxa:
            # get asvs and cluster designation
            asvs = taxdf.loc[taxdf[rank]==tax, "cluster"]
            # generate strings with '<asv> <cluster>' format to match fasta headers
            asvs = [f"{x[0]} {x[1]}" for x in list(zip(asvs.index, asvs))]
            # skip taxa with fewer than 2 ASVs
            if len(asvs) < 2:
                continue
            # write asvs to file
            with open(f"{tax}.ids", "w") as f:
                _ = f.write("\n".join(asvs))
            # use seqkit to extract sequences
            with open(f"{outdir}/{tax}.fasta", 'w') as fhout:
                cmd = ["seqkit", "grep", "-f",f"{tax}.ids","-n", "--quiet", input.fasta]
                _ = subprocess.run(cmd, stdout=fhout)


rule matchlist_vsearch:
    """
    Create a matchlist for sequences in an order using vsearch
    """
    output:
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/vsearch/{tax}.matchlist.tsv"
    input:
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/fasta/{tax}.fasta"
    log:
        "logs/matchlist_vsearch/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{noise_rank}/{tax}.log"
    conda:
        "../envs/vsearch.yml"
    params:
        tmpdir = "$TMPDIR/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/{noise_rank}/{tax}/",
        maxhits = config["noise_filtering"]["max_target_seqs"]
    threads: 1
    shell:
        """
        mkdir -p {params.tmpdir}
        cat {input} | sed 's/>.\+ />/g' > {params.tmpdir}/cluster_reps.fasta
        vsearch --usearch_global {params.tmpdir}/cluster_reps.fasta --db {params.tmpdir}/cluster_reps.fasta --self --id .84 --iddef 1 \
            --userout {output} -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits {params.maxhits} --threads {threads} > {log} 2>&1
        rm -rf {params.tmpdir}
        """

rule generate_datasets:
    """
    Output taxonomy and counts files for a certain taxon
    """
    output:
        taxonomy="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/data/{tax}_taxonomy.tsv",
        counts="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/data/{tax}_counts.tsv",
        cal_counts="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/data/{tax}_cal_counts.tsv",
        tot_prop_counts="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/data/{tax}_tot_prop_counts.tsv",
        sample_prop_counts="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/data/{tax}_sample_prop_counts.tsv",
    input:
        counts=rules.generate_counts_files.output.cluster_counts,
        cal_counts=rules.generate_counts_files.output.calibrated_counts,
        tot_prop_counts=rules.generate_counts_files.output.tot_proportional_counts,
        sample_prop_counts=rules.generate_counts_files.output.sample_proportional_counts,
        taxonomy="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_taxonomy.tsv",
    log:
        "logs/generate_datasets/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{noise_rank}/{tax}.log"
    params:
        split_rank = config["noise_filtering"]["split_rank"]
    conda:
        "../envs/r-env.yml"
    script:
        "../scripts/neeat/generate_datasets.R"

rule generate_aa_seqs:
    """
    Translate nucleotide sequences to amino acid sequences
    """
    output:
        faa="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/faa/{tax}.faa"
    input:
        fasta="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/fasta/{tax}.fasta"
    log:
        "logs/generate_aa_seqs/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{noise_rank}/{tax}.log"
    params:
        codon_table = config["codon_table"],
        codon_start = config["codon_start"]
    conda:
        "../envs/r-env.yml"
    script:
        "../scripts/neeat/generate_aa_seqs.R"

rule mafft_align:
    """
    Align protein sequences using MAFFT
    """
    output:
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/mafft/{tax}.aligned.faa"
    input:
        rules.generate_aa_seqs.output.faa,
    log:
        "logs/mafft_align/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{noise_rank}/{tax}.log",
    conda: 
        "../envs/mafft.yml"
    threads: 4
    shell:
        "mafft --auto --thread {threads} {input} > {output} 2>{log}"

rule pal2nal:
    """
    Generate the corresponding nucleotide alignments with pal2nal
    """
    output:
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/pal2nal/{tax}.fasta"
    input:
        pep=rules.mafft_align.output[0],
        nuc="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/fasta/{tax}.fasta"
    log:
        "logs/noise_filtering/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{noise_rank}/{tax}.log"
    conda:
        "../envs/pal2nal.yml"
    params:
        codon_table = config["codon_table"],
        codon_start = config["codon_start"],
        tmpdir="$TMPDIR/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{noise_rank}/{tax}"
    threads: 1
    shell:
        """
        mkdir -p {params.tmpdir}
        seqkit subseq --region {params.codon_start}:-1 {input.nuc} > {params.tmpdir}/nuc.fasta
        pal2nal.pl {input.pep} {params.tmpdir}/nuc.fasta -output fasta -codontable {params.codon_table} > {output} 2>{log}
        rm -r {params.tmpdir}
        """

rule generate_evodistlists:
    """
    generates order-level evolutionary distance files for the evo_filter function of neeat
    """
    output:
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/evodist/{tax}_evodistlist.tsv"
    input:
        matchlist=rules.matchlist_vsearch.output[0],
        taxonomy=rules.generate_datasets.output.taxonomy,
        fasta=rules.pal2nal.output[0]
    log:
        "logs/generate_evodistlists/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{noise_rank}/{tax}.log"
    params:
        codon_model = workflow.source_path("../scripts/neeat/codon_model.R"),
    conda:
        "../envs/r-env.yml"
    threads: 1
    script:
        "../scripts/neeat/generate_evo_dists.R"

def aggregate_evodistlists(wildcards):
    """
    Aggregate the evodistlists for all taxa
    """
    checkpoint_dir = checkpoints.generate_taxa_seqs.get(**wildcards).output[0]
    files = expand("results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/evodist/{tax}_evodistlist.tsv",
        tool=wildcards.tool,
        rundir=wildcards.rundir,
        chimera_run=wildcards.chimera_run,
        chimdir=wildcards.chimdir,
        rank=wildcards.rank,
        run_name=wildcards.run_name,
        noise_rank=wildcards.noise_rank,
        tax=glob_wildcards(os.path.join(checkpoint_dir, "{tax}.fasta")).tax
        )
    return files

rule evodist:
    output:
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/evodistlist.tsv"
    input:
        aggregate_evodistlists
    shell:
        """
        cat {input} > {output}
        """

rule neeat_filter:
    output:
        "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/filtered/{tax}_filtered.tsv"
    input:
        counts=rules.generate_datasets.output.counts,
        distlist=rules.generate_evodistlists.output[0],
        taxonomy=rules.generate_datasets.output.taxonomy,
    log:
        "logs/neeat_filter/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{noise_rank}/{tax}.log"
    params:
        echo_filter=workflow.source_path("../scripts/neeat/echo_filter.R"),
        evo_filter=workflow.source_path("../scripts/neeat/evo_filter.R"),
        abundance_filter=workflow.source_path("../scripts/neeat/abundance_filter.R"),
    