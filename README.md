[![Pixi Badge](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/prefix-dev/pixi/main/assets/badge/v0.json)](https://pixi.sh)

# HAPP: High-Accuracy Pipeline for Processing deep metabarcoding data

- [Overview](#overview)
- [Installation](#installation)
  - [Software requirements](#software-requirements)
- [How to run the workflow](#how-to-run-the-workflow)
  - [Testrun](#testrun)
  - [Configuration file](#configuration-file)
  - [Configuration profile](#configuration-profile)
  - [Running parts of the workflow](#running-parts-of-the-workflow)
- [Configuration](#configuration)
  - [Taxonomic assignments](#taxonomic-assignments)
  - [Preprocessing](#preprocessing)
  - [Chimera filtering](#chimera-filtering)
  - [ASV clustering tools](#asv-clustering-tools)
  - [Noise filtering](#noise-filtering)
- [Workflow output](#workflow-output)

## Overview

This repository contains a Snakemake workflow for taxonomic assignments, filtering and clustering of Amplicon
Sequence Variants (ASVs). 

Currently the following clustering tools are supported:

| Software  | Reference                                                                                      | Code                                                                    |
|-----------|------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------|
| SWARM     | [Mahé et al 2014](https://peerj.com/articles/593/)                                             | [GitHub](https://github.com/torognes/swarm)                             |
| OptiClust | [Westcott & Schloss 2017](https://journals.asm.org/doi/10.1128/mSphereDirect.00073-17)         | [GitHub](https://github.com/SchlossLab/Westcott_OptiClust_mSphere_2017) |
| dbOTU3    | [Olesen et al 2017](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0176335) | [GitHub](https://github.com/swo/dbotu3)                                 |

The idea with this workflow is to make it easy to run OTU clustering with many 
different parameter settings then evaluate which settings you think works best
for your data.

## Installation

To install the workflow, either clone the repository by running:

```bash
git clone git@github.com:insect-biome-atlas/happ.git
cd happ
```

or visit the [release page](https://github.com/insect-biome-atlas/happ/releases)
and download the latest release.

### Software requirements

The software required to run the workflow can either be installed via
[pixi](https://pixi.sh) or with [Conda](https://docs.conda.io/en/latest/).

#### Installation with pixi (recommended)

Follow the [instructions](https://pixi.sh/latest/#installation) to install pixi
on your system, then run:

```bash
pixi shell
```

from within the root of the repository. This will activate an interactive shell
ready to use with the workflow.

#### Install with Conda

If you prefer to use Conda, you can create a new environment with the required
software by running:

```bash
conda env create -f environment.yml
```

Then activate the environment with:
  
  ```bash
  conda activate happ
  ```

## How to run the workflow

Once you have activated the software environment (either with `pixi shell` or
`conda activate happ` as described above) the basic syntax to run the workflow
is:

```bash
snakemake --sdm <apptainer/conda> \
  --configfile <path-to-your-configfile.yml> \
  --profile <slurm/dardel/local> \
  <additional-arguments>
```

Below is a description of the different command line flags:

- `--sdm <apptainer/conda>`

The `--sdm` flag (short-hand for `--software-deployment-method`) specifies if
rule-specific software dependencies should be handled with
[Apptainer](https://apptainer.org/) (`--sdm apptainer`) or with
[Conda](https://docs.conda.io/en/latest/) (`--sdm conda`). 

> [!NOTE]
> We recommend to use Apptainer if it is available on your system.

- `--configfile <path-to-your-configfile.yml>`

The `--configfile` flag specifies the path to a [configuration
file](#configuration-file) in YAML format. The workflow is preconfigured with
default settings for most parameters (see the `config/config.yml` file for
default settings) but there are a few parameters that you will have to modify in
order to run HAPP on your data. Read more on this in the
[Configuration](#configuration) section.

> [!TIP] 
> We recommend that you make a copy of the `config/config.yml` file,
> modify the copy to fit your data and then supply this file when running HAPP.

- `--profile <slurm/dardel/local/test>`

The `--profile` flag specifies the [configuration
profile](#configuration-profile) to use. These profiles modify the behaviour of
Snakemake itself (such as compute resources etc), in contrast to configuration
files (specified with `--configfile`) which set parameters for HAPP (such as
input files etc).

> [!TIP]
> If you are running the workflow on a high performance computing system with
> the SLURM workload manager use the configuration profile under
> `profiles/slurm`. If you're running on the Dardel HPC system at PDC
> specifically, use the configuration profile under `profiles/dardel`. See the
> README files in the respective subdirectory for more information.

- `<additional-arguments>`

You can append additional Snakemake arguments to the command line call when running HAPP if you wish. These can include number of cores to use (_e.g._ `--cores 4` to run with 4 cores) or performing a dry-run with `-n`. Read about all available Snakemake command line arguments [here](https://snakemake.readthedocs.io/en/stable/executing/cli.html#).

### Testrun

The workflow comes with a relatively small test dataset of 100 sequences in 100
samples which can be used to try the workflow on your system. To use it with
default settings you need to obtain a Sintax reference file and set the
corresponding path in the `sintax:` configuration entry. For ease of use you can
run the `test/setup-ref.sh` script which downloads a reference from
[Figshare](https://doi.org/10.17044/scilifelab.20514192) and updates the
`test/configfile.yml` configuration file:

```bash
bash test/setup-ref.sh
```

Once the above command completes you can run:

```bash
snakemake --profile test --sdm conda
```

to run the workflow with software dependencies handled by Conda, or:

```bash
snakemake --profile test --sdm apptainer
```

to use Apptainer to run jobs in containers where needed.

### Configuration file

The workflow can be configured using a configuration file in YAML format. A
default config file is available at [config/config.yml](config/config.yml). It's
recommended that you copy this file and make your changes in the copy, then pass
it on the command line with `--configfile <path-to-your-configfile>`.

See the [Configuration](#configuration) section for more information on how to configure the workflow.

### Configuration profile

The workflow comes with three different configuration profiles. These are
basically subdirectories containing a `config.yaml` file which sets various
Snakemake command line arguments and define resources for running the workflow
on different systems. The available profiles are:

- `local`: For running the workflow on a local machine (_e.g._ your laptop)
- `slurm`: For running the workflow on a system with the SLURM workload manager
- `dardel`: For running the workflow on the [Dardel HPC system](https://www.pdc.kth.se/hpc-services/computing-systems/dardel-1.1043529) at PDC. 

> [!IMPORTANT]
> You will have to modify the `slurm_account:` setting in the
> `dardel/config.yaml` or `slurm/config.yaml` file in order to use your compute
> account when running on a cluster. For the generic `slurm` profile you will
> also have to change the `slurm_partition:` to match the default partition on
> your cluster.

To use a profile, simply add `--profile <profile-name>` to the Snakemake command line call. For example, to run the workflow on a system with SLURM you would run:

```bash
snakemake --configfile <path-to-your-configfile.yml> \
  --profile slurm \
  --sdm <conda/apptainer>
```

### Running parts of the workflow

Instead of running the workflow all the way through you can also target specific steps such as taxonomic assignments, chimera filtering _etc._. Currently the supported steps are shown in the table below. The general syntax to run up to a certain step is:

```bash
snakemake --configfile <your-configfile.yml> \
  --profile <local/dardel/slurm> \
  --sdm=<conda/apptainer> <step>
```

So to only run up to and including taxonomic assignments for your input sequences, using the default `config/config.yml` configuration file, the `local/` profile and Conda to handle software dependencies you would run:

```bash
snakemake --configfile config/config.yml \
  --profile local \
  --sdm=conda assign_taxonomy
```

> [!IMPORTANT]
> Note the `=` sign in `--sdm=conda`. This is important here in order for Snakemake to properly parse the command line arguments when specific targets are passed.

| Step | Snakemake target | Description |
| ---- | ---------------- | ------- |
| Preprocess | `preprocess` | Performs preprocessing of input sequences as configured under the `preprocessing` section in the config file. |
| Assign taxonomy | `assign_taxonomy` | Assigns taxonomy to the input sequences using tools defined by `taxtools` parameter.
| Filter chimeras | `filter_chimeras` | Filters chimeras using settings defined under `chimera` section in config file. |
| dbotu3 clustering | `dbotu3` | Runs clustering of input sequences using dbotu3. |
| opticlust clustering | `opticlust` | Runs clustering of input sequences using opticlust. |
| swarm clustering | `swarm` | Runs clustering of input sequences using swarm. |
| clustering | `clustering` | Runs clustering of input sequences using all tools defined by `software` parameter in config file. |

### Running the NEEAT algorithm directly

One of the final parts of HAPP is a new noise filtering algorithm called NEEAT
(**N**oise reduction using **E**chos, **E**volutionary signals and **A**bundance
**T**hresholds) which runs on the sequence clusters generated by the workflow. If you want to bypass the other steps of HAPP and run the NEEAT algorithm directly on your sequences you can do so by supplying three required files:

1. A FASTA file containing your sequences. Because HAPP expects the NEEAT algorithm to run on clustered sequences the FASTA headers must be formatted in the pattern `>SEQUENCE_ID CLUSTER_ID` where `SEQUENCE_ID` is the original name of the sequence and `CLUSTER_ID` is the name of the cluster for which the sequence is a representative. When plugging in your data you can use the same name for both `SEQUENCE_ID` and `CLUSTER_ID`, _e.g._ `>sequence1 sequence1`.
2. A tab-separated file with taxonomic information about each sequence. This file must have the sequence ids in the first column and must contain the columns `cluster` and `representative` where the `cluster` column contains the name of the cluster to which the sequence was assigned and `representative` contains a `1` if the sequence is a representative of the cluster, otherwise it contains a `0`. Again, you may simply set the `cluster` column to the same name as the first column and set the `representative` column to `1` for all your sequences when plugging your data into NEEAT this way. By default the workflow partitions the sequences by taxonomy and runs NEEAT in parallell on each partition for better efficiency. The default rank at which this partition occurs is `Order` and is configured via the `split_rank` parameter under the `noise_filtering` [section](#noise-filtering) in the configuration file. You can modify this parameter but make sure that the value you set for `split_rank` matches with a column in this tab-separated file.
3. A tab-separated file with counts of sequences (rows) in each sample (columns). The names in the first column must match with the `CLUSTER_ID` in the headers of the FASTA file (see point 1), and with the `cluster` column in the tab-separated information file (see point 2).

To see an example of these three files take a look at the `cluster_reps.fasta`, `cluster_taxonomy.tsv` and `cluster_counts.tsv` files in the [data/neeat_test/](https://github.com/insect-biome-atlas/happ/tree/main/data/neeat_test) directory supplied with the repository.

These files are specified as input to the standalone version of NEEAT by adding a `neeat:` section to your configuration file, like so:

```yaml
neeat:
  fastafile: <path-to-point1-fastafile>
  taxfile: <path-to-point2-TSVfile>
  countsfile: <path-to-point3-TSVfile>
```

See the `test/configfile.yml` file supplied with the workflow for an example.

Once you have got your files and added the parameters to the configuration file you can run NEEAT directly on your data with:

```bash
snakemake --configfile <your-configfile.yml> \
  --profile <local/dardel/slurm> \
  --sdm conda \
  -s workflow/rules/neeat-standalone.smk
```

The output will be placed under `results/neeat/` and the main results files will be in a subdirectory with the name of the taxonomic rank used to partition the data. So by default there will be a directory `results/neeat/Order`. This directory contains results from intermediate steps in different subdirectories and the files:

- `discarded_cluster_taxonomy.tsv`: Information about sequences discarded as noise by NEEAT.
- `noise_filtered_cluster_taxonomy.tsv`: Information about sequences retained after NEEAT filtering.
- `noise_filtered_cluster_counts.tsv`: Counts of sequences retained after NEEAT filtering.

To test NEEAT on a small dataset supplied with the workflow you may run:

```bash
snakemake --profile test \
  --sdm conda \
  -s workflow/rules/neeat-standalone.smk
```


## Configuration

The most important parameters in the config file are the input files:

- `asv_seqs.fasta` (ASV sequences in FASTA format)
- `asv_counts.tsv` (tab separated file with counts of ASVs (rows) in samples (columns))

These files **must** be placed in a subdirectory under `data/` that is specified by the `rundir` parameter in the config file. The `rundir` parameter should be the name of the subdirectory containing the input files.

The default config file contains:

```yaml
rundir: "test"
```

and the input files are placed under `data/test/`.

The `run_name:` parameter defines the name of the run and can be used to
separate different runs of the workflow (_e.g._ with different settings on the
same input data). The default is `run1`.

The `split_rank:` parameter specifies the taxonomic rank to split the ASVs by
before running the ASV clustering tools. Splitting the data means that ASVs that
do not share the same `split_rank` are not compared for clustering which means
that you should set this parameter to a relatively high rank in the taxonomy
tree. The default is `Family`.

The `ranks:` parameter 

### Taxonomic assignments

HAPP uses taxonomic information about your input ASVs to split sequences by
taxonomic ranks, allowing for parallel clustering of ASVs. Therefore the ASVs
need to be assigned to taxonomic ranks. This can be done using the workflow
itself, or by providing a file with taxonomic assignments.

The workflow supports taxonomic assignment using the following tools:

- SINTAX
- EPA-NG (including assignments with GAPPA)
- SINTAX + EPA-NG (reassigning SINTAX assignments with EPA-NG)
- vsearch (implemented in QIIME2)
- sklearn (implemented in QIIME2)

To assign taxonomy to your input ASVs you need to provide a reference database
that works with SINTAX and/or with EPA-NG and GAPPA, depending on the tools you
want to use. See below for instructions on how to obtain and configure these
references.

> [!IMPORTANT] 
> You will also need to make sure that the `split_rank` and `ranks`
> config parameters are set correctly to match the taxonomic ranks in your
> taxonomic assignments file.

#### SINTAX

A SINTAX compatible reference database containing COI (mitochondrial cytochrome oxidase subunit I) sequences collected from the BOLD database is available on [Figshare](https://doi.org/10.17044/scilifelab.20514192). From there, simply download the file `bold_clustered.sintax.fasta.gz`, unzip it and set the path to the file in the config file under the `sintax` section:

```yaml
sintax:
  ref: "/path/to/bold_clustered.sintax.fasta"
```

By default, SINTAX assignments are made using a cutoff of 0.8. This can be changed by setting the `cutoff` parameter in the config file under the `sintax` section.

#### EPA-NG and GAPPA assignments

The phylogenetic placement/assignment tools EPA-NG and GAPPA require a reference
tree, a multiple alignment and a reference taxonomy file. A compatible reference that allows assignments to classes Collembola, Diplura, Protura and Insecta is available to download from https://github.com/insect-biome-atlas/paper-bioinformatic-methods/tree/main/data/chesters_tree. The files you need are:

- chesters_new_outgroups_aligned.trim0.9.fasta (alignment)
- chesters_new_outgroups.nwk (tree)
- taxonomy.tsv (taxonomy)

Download these files and place them in a directory on your system. Then set the path to the files in the config file under the `epa-ng` section in the config file:

```yaml
epa-ng:
  tree: "/path/to/chesters_tree/chesters_new_outgroups.nwk"
  msa: "/path/to/chesters_new_outgroups_aligned.trim0.9.fasta"
  ref_taxonomy: "/path/to/chesters_tree/taxonomy.tsv"
```

EPA-NG and GAPPA can be further configured with the following parameters under the `epa-ng` section in the config file:

- `model:` This is the starting model used by RAxML-NG to evaluate the reference tree.
- `heuristic:` This is the preplacement heuristic used by EPA-NG. The default is `dyn-heur` which corresponds to the default dynamic mode in EPA-NG. Also available are `baseball-heur` (as used by pplacer) and `no-heur` (no preplacement heuristic).
- `chunk_size:` This is the number of sequences to process in at a time. Set this to a lower value if you run out of memory during the EPA-NG step.

Gappa can be configured using the `gappa:` section under `epa-ng:`:

- `distribution_ration:` Determines how gappa handles edges with two possible annotations. If set to `-1` (default) the program will determine the ratio automatically from the 'distal length' specified per placement.
- `consensus_thresh:` When assigning taxonomy to missing labels in the reference, require this consensus threshold to assign a label. The default is `1`.

#### VSEARCH/sklearn assignments

To run the `vsearch` or `sklearn` taxonomic assignment methods (implemented in QIIME2) you need to provide paths to files compatible with QIIME2 under the `qiime2:` config section. 

- `ref:` This is the path to a fasta file containing reference sequences for the taxonomic assignment with the QIIME2 `feature-classifier` plugin.
- `taxfile:` This is the path to a file containing the taxonomy for the reference sequences in the fasta file.
- `ranks:` This is a list of taxonomic ranks present in the reference.

See the [QIIME2 docs](https://amplicon-docs.qiime2.org/en/latest/references/plugins/feature-classifier.html) for how to obtain and format these files.

> [!TIP]
>If you downloaded the COIDB SINTAX reference above you can leave the `ref:`,
>`taxfile:` and `ranks:` config parameters empty. HAPP will then generate the
>necessary files using the SINTAX reference.

#### Combining SINTAX + EPA-NG assignments

If you want to run SINTAX and EPA-NG then use the latter to reasssign SINTAX assignments include `sintax+epa-ng` in the `taxtools` list in the config file:

```yaml
taxtools:
  - "sintax+epa-ng"
```

To use this as the taxonomic source for downstream steps, also set the `taxonomy_source:` parameter in the config file to `sintax+epa-ng`:

```yaml
taxonomy_source: "sintax+epa-ng"
```

We have found that SINTAX has a high precision at low taxonomic ranks (_e.g._ species) meaning when it makes an assignment those assignments are often correct. However, SINTAX is conservative meaning that it will often fail to make assignments at higher taxonomic ranks (_e.g._ order), leaving ASVs unclassified. EPA-NG on the other hand is less conservative and will often make correct assignments at higher ranks such as order. By combining the two methods we can take advantage of the strengths of both, using SINTAX to make assignments at low ranks and use EPA-NG to update higher level assignments if those exist.

The `reassign:` config section handles this behaviour and can be configured with the following parameters:

- `placeholder_rank:` This is the taxonomic rank for which both SINTAX and EPA-NG assignments must be identical in order for any updates to be made. The default is `Class` which means that if both SINTAX and EPA-NG agree on the class level assignment for an ASV, the EPA-NG assignment can be used to update the SINTAX assignment (if additional criteria are met, see below).
- `placeholder_taxa:` This is a list of taxa names at rank = `placeholder_rank` for which to make updates. This list should only contain taxa with sufficient coverage in the reference database used for EPA-NG.
- `reassign_ranks:` This is a list of taxonomic ranks at which to update assignments. The default is `["Order"] which means that assignments will be updated at the order level if all criteria are met.
- `downstream_ranks:` This is a list of taxonomic ranks below the `reassign_ranks`. Taxonomic assignments for ranks in this list will be prefixed with 'unclassified.' if the assignment at `reassign_ranks` is updated. The default is `["Family", "Genus", "Species", "BOLD_bin"]`.

With default settings, this means that ASVs without an assignment at Order in SINTAX will be updated with the Order level EPA-NG assignment if the EPA-NG assignment is identical to the SINTAX assignment at Class. Downstream ranks will be left unchanged, unless they are also unassigned in which case they will be prefixed with 'unclassified.' followed by the updated assignment at the rank above.

#### Existing taxonomic assignments (optional)

If you already have taxonomic assignments for your ASVs you can point to the file with taxonomic assignments in the config file with the `taxonomy_source:` parameter:

```yaml
taxonomy_source: "/path/to/asv_taxa.tsv"
```

 This file should be a tab-separated file with the ASV id in the first column and subsequent columns with taxonomic assignments at different ranks. An example of such a file is shown below:

| ASV | Kingdom | Phylum | Class | Order | Family | Genus | Species | BOLD_bin |
|-----|---------|--------|-------|-------|--------|-------|---------|----------|
| ASV1 | Animalia | Arthropoda | Insecta | Coleoptera | Scarabaeidae | Melinopterus |  Melinopterus punctatosulcatus | BOLD:AAJ8574 |
| ASV2 | Animalia | Arthropoda | Insecta | Lepidoptera | Tortricidae | Dichrorampha | Dichrorampha sylvicolana | BOLD:AAL6971 |

If you do not want to run additional taxonomic assignment tools, you can then set the `taxtools:` parameter to an empty list in the config file:

```yaml
taxtools: []
```

#### Consensus taxonomy

Once the ASVs have been clustered into OTU clusters, the workflow will attempt to assign a consensus taxonomy to each cluster. This is done by comparing the taxonomic assignments of the ASVs in the cluster starting from the lowest taxonomic rank and moving up until a consensus is reached. The consensus is reached when the sum of the abundance weighted assignments for a taxon is equal to or greater than a certain threshold (configurable by the `consensus_threshold` parameter in the config, default is 80%). The consensus taxonomy is then assigned to the cluster. The `consensus_ranks` parameter in the config file specifies the taxonomic ranks to use when assigning consensus taxonomy. The default is `["Family", "Genus", "Species"]` which means that the workflow will only attempt to assign a consensus using these ranks.

> [!IMPORTANT]
> The ranks specified with `consensus_ranks:` must be present in the taxonomic
> assignments file used as input to the workflow. If a rank is missing from the
> taxonomic assignments file, the workflow will not be able to assign a
> consensus taxonomy at that rank. If you are using a custom taxonomic
> assignments file, make sure that it contains the ranks specified in the
> `consensus_ranks:` parameter.

As an example, if a cluster contains 3 ASVs with the following taxonomic
assignments and total sum of counts across samples:

| ASV | Family | Genus | Species | ASV_sum |
|-----|--------|-------|---------|---------|
| ASV1 | Tenthredinidae | Pristiphora | Pristiphora mollis | 20 |
| ASV2 | Tenthredinidae	| Pristiphora	| Pristiphora cincta | 60 |
| ASV3 | Tenthredinidae	| Pristiphora	| Pristiphora leucopodia | 20 |

then at Species the abundance weighted taxonomic assignment is 60% _Pristiphora
cincta_ and 20% _Pristiphora mollis_ and _Pristiphora leucopodia_ each. At a 80%
consensus threshold we cannot assign a taxonomy at Species level to the cluster,
so the algorithm will move up to Genus level where we have 100% agreement for
Pristiphora. The final consensus taxonomy for the cluster will then be:

| Cluster | Family | Genus | Species |
|---------|--------|-------|---------|
| Cluster1 | Tenthredinidae | Pristiphora | unresolved.Pristiphora |

At `consensus_threshold: 60` we would have been able to assign taxonomy at the
Species level and the cluster taxonomy would have been:

| Cluster | Family | Genus | Species |
|---------|--------|-------|---------|
| Cluster1 | Tenthredinidae | Pristiphora | Pristiphora cincta |

### Preprocessing

The workflow supports optional preprocessing of the input data by filtering ASVs by length and/or removal of ASVs with in-frame stop codons. These steps are controlled by the `preprocess` section in the config file.

The following parameters are available:

- `filter_length:` A boolean specifying whether to filter ASVs by length. The default is `False`.
- `min_length:` The minimum length of ASVs to keep. The default is `403`.
- `max_length:` The maximum length of ASVs to keep. The default is `418`.
- `filter_codons:` A boolean specifying whether to filter ASVs for in-frame stop codons. The default is `False`.
- `stop_codons:` A comma separated list of stop codons to look for. The default is `"TAA,TAG"`.
- `start_position:` The position in the ASV sequence to start looking for stop codons. The default is `2`.
- `end_position:` The position in the ASV sequence to stop looking for stop codons. The default is `0` which means that the entire sequence is checked.

### Chimera filtering

The workflow supports optional chimera removal using the uchime algorithm 
implemented in vsearch. Chimera detection can be run either in 'batchwise' mode
using the data under `rundir` directly or in 'samplewise' mode in which the 
`asv_counts.tsv` file is used to generate sample-specific fasta files containing 
sequences with a count > 0 in each sample. Sequences identified as non-chimeric
are then passed downstream in the workflow and used as input for clustering.

Config parameters for the chimera removal part of the workflow are nested 
under the `chimera:` section in the config file:

- `remove_chimeras:` A boolean specifying whether to run chimera filtering or not. The default is `True`.
- `run_name:` A name for the chimera filtering run. This is used to separate different runs of the chimera filtering. The default is `chimera1`.
- `method:` The chimera filtering method to use. Can be either `batchwise` or `samplewise`. The default is `samplewise`.
- `algorithm:` The chimera detection algorithm to use. Can be either
`uchime_denovo`, `uchime2_denovo` or `uchime3_denovo`. The default is
`uchime_denovo`. The latter two require perfect matches between the ASV sequence
and a chimeric model, whereas 'uchime_denovo' does not.
- `min_samples_shared:` The number of samples in which chimeric ASVs have to be present with their 'parents' in order to be marked as chimeric. The default is `1`.
- `min_samples_shared:` parameter is specific to the batchwise mode and
specifies the number of samples in which chimeric ASVs have to be present with
their 'parents' (see Uchime
[docs](https://drive5.com/usearch/manual/chimeras.html)) in ordered to be marked
as chimeric.
- `min_frac_samples_shared:` parameter is similar to `min_samples_shared`, but
instead of an absolute number require that sequences are present with their
parents in a fraction of the samples in which they are present.
- `min_chimeric_samples:` refers to the samplewise chimera mode, and
determines in how many samples an ASVs has to be marked as chimeras in to be
filtered as chimeric. Setting this value to 0 (the default) means that sequences
have to be marked as chimeric in **all** samples where they are present in order
to be filtered as chimeric.
- `min_frac_chimeric_samples:` is analogous to `min_chimeric_samples` but specifies
a fraction of samples, instead of an absolute number

The `dn`, `mindiffs`, `mindiv` and `minh` parameters are specific to how the
chimeric score of sequences is calculated. Please see the Uchime
[docs](https://www.drive5.com/usearch/manual6/UCHIME_score.html) for details. 

In all algorithms ASV sequences are first aligned to other ASVs that are above a
certain abundance threshold. This so called [abundance
skew](https://www.drive5.com/usearch/manual6/abundance_skew.html) threshold is
by default set to 2.0 for the 'uchime_denovo' and 'uchime2_denovo' algorithms
and to 16.0 for 'uchime3_denovo'. However, this value can be overridden by
explicitly setting `abskew` in the config file under `chimera:`.

### ASV clustering tools

The workflow can be configured to run one or more tools to cluster ASVs into
OTUs. The `software` parameter takes a list of supported tool names: `swarm`,
`opticlust` and `dbotu3` currently.

Below are specific parameters for tools used in the workflow.

#### vsearch

Pairwise identities between ASVs are obtained by running vsearch with the
`usearch_global` setting. This is used as input by `opticlust` for clustering
ASVs. Below are the parameters specific to vsearch:

- `id:` sets the minimum pairwise identity to report. The default is 0.84 meaning that sequences with a pairwise identity less than 84% is rejected.
- `iddef:` sets the identity definition used by vsearch. The default
setting of `1` means that identity is calculated as `(matching columns) /(alignment length)`. Refer to the vsearch manual for other settings.
- `query_cov:` sets the minimum aligned fraction of the query sequence. The
default setting of `0.9` means that alignments are rejected if this fraction is
less than 90%.

#### Opticlust

The opticlust method is implemented in `mothur`. Here the optimal clustering of 
sequences is found by iteratively moving sequences between OTU clusters in an 
attempt to maximize the Matthew’s Correlation Coefficient (MCC) which is a metric
that combines True/False positives and negatives. In this context:
- true positives: are sequences that are within a maximum pairwise distance from each
other (defined by the `cutoff` parameter) and that share an OTU cluster
- false positives are further in distance than the cutoff and share a cluster,
- false negatives are within the cutoff but are placed in different clusters, and
- true negatives are further in distance than the cutoff and are placed in different clusters

Opticlust takes as input a pairwise similarity matrix and a file with sum of counts
across samples for each ASV.

For a full description of opticlust parameters, see the
[manual](https://mothur.org/wiki/cluster/). Below are the parameters specific to
opticlust:

- `aligner:` specifies whether to use `vsearch` (default) or the built-in
aligner in `mothur` for generating pairwise identities. Note that the latter has
not been tested for this workflow.
- `delta:` sets the stable value for the metric in opticlust.
- `cutoff:` controls at what similarity threshold clusters are generated. 
- `initialize:` sets the initial randomization for opticlust. The default is
`singleton` where each sequence is randomly assigned to its own OTU. The other
accepted setting `oneotu` means that all sequences are assigned to one otu. 
- `precision:` sets the floating point precision for opticlust.

#### Swarm

Swarm clusters sequences using a local linking threshold `d` which represents
the maximum number of differences between two amplicons. The input is a single
fasta sequence where the total count of each sequence is suffixed to the
sequence ID.

For a full description of Swarm parameters, see the
[manual](https://github.com/torognes/swarm/blob/master/man/swarm_manual.pdf).
Below are the parameters specific to Swarm:

- `differences:` specifies the maximum number of differences (mismatches)
allowed between two ASVs, meaning that two ASVs will be grouped if they have
integer (or less) differences.
- `no-otu-breaking:` deactivates Swarms built-in cluster refinement
when using `d=1`.
- `fastidious:` specifies whether to perform a second clustering pass
(when set to True) to reduce the number of small clusters when using `d=1`.
- `boundary:` defines the minimum abundance of what should be considered a
_large_ cluster when using `d=1`.

The parameters `match-reward`, `mismatch-penalty`, `gap-opening-penalty` and
`gap-extension-penalty` allows control of Swarms advanced Pairwise alignment
options.

#### dbotu3

dbotu3 uses both sequences and their distribution across samples to cluster ASVs
into OTUs. It takes as input ASV sequences in fasta format and a table of
counts.

For a full description of dbotu3 parameters, see the
[documentation](https://dbotu3.readthedocs.io/en/latest/). Below are the
parameters specific to dbotu3:

- `dist:` sets the maximum allowed genetic dissimilarity between
sequences.
- `abund:` sets the minimum fold difference for comparing two OTUs.
- `pval:` sets the minimum p-value for merging ASVs.

### Noise filtering

The generated clusters can be filtered for noise such as Nuclear mitochondrial
DNA segments (NUMTs) and sequencing errors. This is done using taxonomic
information, sequence similarity and abundance information. Below are the
parameters specific to the noise filtering step:

- `run_name:` sets the name of the run for the noise filtering step. Similarly
  to the `chimera` section, this is used to separate different runs of the noise
  filtering.
- `split_rank:` sets the taxonomic rank to split the ASVs by before running the
  noise filtering.
- `assignment_rank:` sets the taxonomic rank to use for the noise filtering.
  ASVs unassigned or ambiguously assigned at this rank are considered noise and
  will be removed. To skip this step set `assignment_rank: ""`.
- `max_target_seqs:` sets the maximum number of target sequences to return for the vsearch alignent step.
- `min_match:` sets the minimum threshold of sequence similarity (in %) for
  considering any OTU as an error of another.
- `n_closest:` sets the number of potential parent clusters to use in the
  comparison when identifying noise.
- `echo_min_overlap:` in the 'echo filter', sets the minimum fraction of the
  putative “noise” OTU samples that must also contain the putative authentic CO1
  “parent” OTU.
- `echo_max_read_ratio:` sets the maximum value of the ratio between noise and
  parent reads.
- `echo_read_ratio_type:` sets the type of read ratio considered by the
  `echo_max_read_ratio` criterion. Can be either `max` or `mean`.
- `echo_require_corr:` specifies whether a significant correlation between noise
  and parent read numbers should be required before the putative OTU is removed.
- `evo_local_min_overlap:` in the 'local evo filter', sets the minimum fraction
  of the putative “noise” OTU samples that must also contain the putative
  authentic CO1 “parent” OTU.
- `dist_type_local:` in the 'local evo filter', specifies whether distances 
  between representative ASVs should be based on unweighted amino acid distances
  (`dadn`) or amino acid distances weighted based on biochemical similarities
  (`wdadn`).
- `dist_threshold_local:` in the 'local evo filter', sets the distance value
  above which an OTU is considered to represent noise.
- `dist_threshold_global:` in the 'global evo filter', sets the distance value
  above which an OTU is considered to represent noise.
- `abundance_cutoff_type:` specifies whether the abundance cutoff should be based
  on the sum of reads across samples in the dataset (`sum`), or the max number of
  reads (`max`).
- `abundance_cutoff:` sets the minimum number of reads required for an OTU to
  pass the abundance filter.
- `codon_table:` sets the codon table to use for translation when aligning
  sequences with pal2nal.
- `codon_start:` sets the reading frame start (1 based).

## Workflow output

HAPP generates results in a directory `results/` with sub-directories organised by several of the parameters set in the configuration used. This allows you to run HAPP with several configurations in the same root directory without overwriting results.

### Preprocessing

The output from preprocessing, if configured to run in your configuration file, is found in `results/preprocess/<rundir>`. If using the test dataset (`data/test/`) and configuration file (`test/configfile.yml`) the output will be:

```
results/preprocess
└── test # <- rundir parameter
    ├── ASV_codon_filtered.fna # <- fasta file after filtering by stop codons
    ├── ASV_codon_filtered.list # <- list of sequences removed by stop codon filtering
    ├── ASV_codon_filtered.table.tsv # <- counts file after filtering by stop codons
    ├── ASV_length_filtered.fna # <- fasta file after filtering by length
    └── ASV_length_filtered.table.tsv # <- counts file after filtering by length
```

If both `filter_length` and `filter_codons` are set to `True`, length filtering happens before stop codon filtering.

### Taxonomic assignments

The output from taxonomic assignments is found in `results/taxonomy` with
sub-directories for each tool used, _e.g._:

```
results/taxonomy
├── epa-ng # output from epa-ng tool
│   └── test # name of rundir
│       └── assignments
│           └── dyn-heur # epa-ng placement heuristic
│               └── taxonomy.tsv # taxonomic assignments
├── sintax # output from sintax tool
│   └── test # name of rundir
│       ├── confidence.tsv # TSV file with sintax confidence values per sequence/rank
│       └── taxonomy.tsv # taxonomic assignments
├── sintax_epang # output from combining sintax + epa-ng
│   └── test # name of rundir
│       └── dyn-heur # epa-ng placement heuristic
│           └── taxonomy.tsv # taxonomic assignments
└── vsearch # output from vsearch tool
    └── test # name of rundir
        ├── taxonomy.raw.tsv # raw output from vsearch, incl. confidence values
        └── taxonomy.tsv # taxonomic assignments
```

### Chimera filtering

If chimera filtering is enabled results from filtering are placed under `results/chimera/<rundir>` with additional sub-directories depending on the `run_name`, `method` and `algorithm` parameter settings under the `chimera` section in the configuration file. For example if running chimera filtering with the `samplewise` method (default):

```
results/chimera/
└── test # <- rundir
    └── filtered
        └── chimera1 # <- run_name under chimera config section
            └── samplewise.uchime_denovo # <- chimera method.algorithm settings
                ├── chimeras.fasta # <- fasta file with chimeric sequences
                └── nonchimeras.fasta # <- fasta file with non-chimeric sequences
```

If chimera filtering is run with the `samplewise` method then each sample will also have a sub-directory under `results/chimera/<rundir>/samplewise.<algorithm>/samples/<sample>` with intermediate files from the chimera filtering steps. For example:

```
results/chimera/
└── test # name of rundir
    └── samplewise.uchime_denovo # method.algorithm chimera settings
        └── samples
            └── sample1 # output for sample1
                ├── borderline.fasta.gz # sequences marked as 'borderline' chimeras
                ├── chimeras.fasta.gz # chimeric sequences
                ├── nonchimeras.fasta.gz # non-chimeric sequences
                ├── uchimealns.out.gz # alignment output file
                └── uchimeout.txt.gz # chimera results file
```

If instead the chimera filtering is run with the `batchwise` method, there will be:

```
results/chimera/
└── test
    ├── batchwise.uchime_denovo
    │   ├── borderline.fasta
    │   ├── chimeras.fasta
    │   ├── nonchimeras.fasta
    │   ├── uchimealns.out
    │   └── uchimeout.txt
    └── filtered
        └── chimera1
            └── batchwise.uchime_denovo
                ├── chimeras.fasta
                ├── nonchimeras.fasta
                └── uchimeout.tsv
```

### Output from clustering tools

Each clustering tool specified by the `software` config parameter gets a separate sub-directory under `results/`, _e.g._ `results/swarm`, which is further structured by the following parameters:

```yaml
rundir: <name of rundir>
split_rank: <taxonomic rank at which to partition input>
run_name: <main name of the workflow run>
chimera:
  run_name: <name of chimera run>
  method: <samplewise/batchwise>
  algorithm: <uchime_denovo/uchime_denovo2/uchime_denovo3>
```

For example, using the `test/configfile.yml` configuration file the results from clustering with `swarm` would be:

```
results/swarm
└── test # rundir
    └── chimera1 # chimera run_name
        └── samplewise.uchime_denovo # chimera method.algorithm
            └── Family # split_rank
                └── runs
                    └── testrun # run_name
                        ├── cluster_consensus_taxonomy.tsv # consensus taxonomy of clusters
                        ├── cluster_counts.tsv # summed counts of clusters
                        ├── cluster_reps.fasta # representative sequences for clusters
                        ├── cluster_taxonomy.tsv # taxonomic info and cluster ass
                        ├── neeat # output from noise filtering with NEEAT
                        ├── precision_recall.order.txt # precision/recall values per Order
                        └── precision_recall.txt # precision/recall values for clustering results
```