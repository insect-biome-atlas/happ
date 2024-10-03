localrules:
    generate_taxa_seqs,
    taxonomy_filter,
    generate_aa_seqs,

rule taxonomy_filter:
    output:
        taxonomy="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/taxonomy_filter/{assignment_rank}/cluster_taxonomy.tsv",
    input:
        taxonomy="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_taxonomy.tsv"
    log:
        "logs/taxonomy_filter/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{assignment_rank}.log"
    run:
        import pandas as pd
        assignment_rank = wildcards.assignment_rank
        taxdf = pd.read_csv(input.taxonomy, sep="\t", index_col=0)
        taxdf = taxdf.loc[(~taxdf[assignment_rank].str.contains("_X+$"))&(~taxdf[assignment_rank].str.startswith("unclassified"))]
        taxdf.to_csv(output.taxonomy, sep="\t")

def get_taxonomy(wildcards):
    if config["noise_filtering"]["assignment_rank"] == "":
        return expand("results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_taxonomy.tsv",
            tool=wildcards.tool, rundir=wildcards.rundir, chimera_run=wildcards.chimera_run, chimdir=wildcards.chimdir, rank=wildcards.rank, run_name=wildcards.run_name)
    else:
        return expand("results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/taxonomy_filter/{assignment_rank}/cluster_taxonomy.tsv",
            tool=wildcards.tool, rundir=wildcards.rundir, chimera_run=wildcards.chimera_run, chimdir=wildcards.chimdir, rank=wildcards.rank, run_name=wildcards.run_name, assignment_rank=config["noise_filtering"]["assignment_rank"])

rule generate_counts_files:
    """
    Generate the counts files for use with neeat filtering
    """
    output:
        cluster_counts="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/counts/cluster_counts.tsv",
        #calibrated_counts="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/counts/calibrated_cluster_counts.tsv",
        #tot_proportional_counts="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/counts/tot_proportional_cluster_counts.tsv",
        #sample_proportional_counts="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/counts/sample_proportional_cluster_counts.tsv",
    input:
        counts="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_counts.tsv",
        taxonomy=get_taxonomy,
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
        taxonomy=get_taxonomy,
        fasta="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/cluster_reps.fasta",
        counts=rules.generate_counts_files.output.cluster_counts
    log:
        "logs/generate_taxa_seqs/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{noise_rank}.log"
    shadow: "minimal"
    run:
        import subprocess
        import pandas as pd
        rank = wildcards.noise_rank
        outdir=output[0]
        os.makedirs(outdir, exist_ok=True)
        # read in the taxonomy file
        taxdf = pd.read_csv(input.taxonomy[0], sep="\t", index_col=0)
        # extract representative ASVs
        taxdf = taxdf.loc[taxdf["representative"]==1]
        # also make sure that only clusters that remain in the counts file are used
        clusters_w_counts = []
        with open(input.counts, 'r') as fhin:
            for line in fhin:
                clusters_w_counts.append(line.split("\t")[0])
        taxdf = taxdf.loc[taxdf["cluster"].isin(clusters_w_counts)]
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
        #cal_counts="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/data/{tax}_cal_counts.tsv",
        #tot_prop_counts="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/data/{tax}_tot_prop_counts.tsv",
        #sample_prop_counts="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/data/{tax}_sample_prop_counts.tsv",
    input:
        counts=rules.generate_counts_files.output.cluster_counts,
        #cal_counts=rules.generate_counts_files.output.calibrated_counts,
        #tot_prop_counts=rules.generate_counts_files.output.tot_proportional_counts,
        #sample_prop_counts=rules.generate_counts_files.output.sample_proportional_counts,
        taxonomy=get_taxonomy,
    log:
        "logs/generate_datasets/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{noise_rank}/{tax}.log"
    params:
        split_rank = config["noise_filtering"]["split_rank"]
    conda:
        "../envs/r-env.yml"
    threads: 1
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
    generates evolutionary distance files for the evo_filter function of neeat
    """
    output:
        tsv="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/evodist/{tax}_evodistlist.tsv"
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

def agg_evodist(wc):
    """
    Aggregate the evodistlists for all taxa
    """
    checkpoint_dir = checkpoints.generate_taxa_seqs.get(**wc).output[0]
    return expand("results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/evodist/{tax}_evodistlist.tsv",
        tool=wc.tool, rundir=wc.rundir, chimera_run=wc.chimera_run, chimdir=wc.chimdir, rank=wc.rank,
        run_name=wc.run_name, noise_rank=wc.noise_rank, tax=glob_wildcards(os.path.join(checkpoint_dir, "{tax}.fasta")).tax)

rule aggregate_evodist:
    input:
        agg_evodist
    output:
        touch("results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/evodist.done")

rule generate_neeat_filtered:
    output:
        retained="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/filtered/{tax}_retained.tsv",
        discarded="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/filtered/{tax}_discarded.tsv",
        counts="results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/filtered/{tax}_counts.tsv",
    input:
        counts=rules.generate_datasets.output.counts,
        distlist=rules.generate_evodistlists.output[0],
        taxonomy=rules.generate_datasets.output.taxonomy,
    log:
        "logs/neeat_filter/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/{run_name}/{noise_rank}/{tax}.log"
    params:
        neeat_filter=workflow.source_path("../scripts/neeat/neeat_filter.R"),
        echo_filter=workflow.source_path("../scripts/neeat/echo_filter.R"),
        evo_filter=workflow.source_path("../scripts/neeat/evo_filter.R"),
        abundance_filter=workflow.source_path("../scripts/neeat/abundance_filter.R"),
        min_match=config["noise_filtering"]["min_match"],
        n_closest=config["noise_filtering"]["n_closest"],
        echo_min_overlap=config["noise_filtering"]["echo_min_overlap"],
        echo_read_ratio_type=config["noise_filtering"]["echo_read_ratio_type"],
        echo_max_read_ratio=config["noise_filtering"]["echo_max_read_ratio"],
        echo_require_corr=config["noise_filtering"]["echo_require_corr"],
        evo_local_min_overlap=config["noise_filtering"]["evo_local_min_overlap"],
        dist_type_local=config["noise_filtering"]["dist_type_local"],
        dist_threshold_local=config["noise_filtering"]["dist_threshold_local"],
        dist_threshold_global=config["noise_filtering"]["dist_threshold_global"],
        abundance_cutoff_type=config["noise_filtering"]["abundance_cutoff_type"],
        abundance_cutoff=config["noise_filtering"]["abundance_cutoff"],
        assignment_rank=config["noise_filtering"]["assignment_rank"]
    conda:
        "../envs/r-env.yml"
    script:
        "../scripts/neeat/generate_neeat_filtered.R"
    
def aggregate_neeat(wc):
    """
    Aggregate the evodistlists for all taxa
    """
    checkpoint_dir = checkpoints.generate_taxa_seqs.get(**wc).output[0]
    retained = expand("results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/filtered/{tax}_retained.tsv",
        tool=wc.tool, rundir=wc.rundir, chimera_run=wc.chimera_run, chimdir=wc.chimdir, rank=wc.rank,
        run_name=wc.run_name, noise_rank=wc.noise_rank, tax=glob_wildcards(os.path.join(checkpoint_dir, "{tax}.fasta")).tax)
    discarded = expand("results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/filtered/{tax}_discarded.tsv",
        tool=wc.tool, rundir=wc.rundir, chimera_run=wc.chimera_run, chimdir=wc.chimdir, rank=wc.rank,
        run_name=wc.run_name, noise_rank=wc.noise_rank, tax=glob_wildcards(os.path.join(checkpoint_dir, "{tax}.fasta")).tax)
    counts = expand("results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/filtered/{tax}_counts.tsv",
        tool=wc.tool, rundir=wc.rundir, chimera_run=wc.chimera_run, chimdir=wc.chimdir, rank=wc.rank,
        run_name=wc.run_name, noise_rank=wc.noise_rank, tax=glob_wildcards(os.path.join(checkpoint_dir, "{tax}.fasta")).tax)
    return {"retained": retained, "discarded": discarded, "counts": counts}

rule neeat:
    output:
    #    counts = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/filtered_counts.tsv",
    #    retained = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/cluster_taxonomy_retained.tsv",
    #    discarded = "results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/cluster_taxonomy_discarded.tsv"
        touch("results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/neeat/{noise_rank}/neeat.done")
    input:
        unpack(aggregate_neeat)
    #run:
    #    import pandas as pd
    #    counts = pd.DataFrame()
    #    for f in input.counts:
    #        _counts = pd.read_csv(f, sep="\t", index_col=0)
    #        _counts = _counts.loc[:, sorted(counts.columns)]
    #        counts = pd.concat([counts, _counts], axis=0)