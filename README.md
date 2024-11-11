[![Pixi Badge](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/prefix-dev/pixi/main/assets/badge/v0.json)](https://pixi.sh)

# HAPP: High-Accuracy Pipeline for Processing deep metabarcoding data

- [Overview](#overview)
- [Installation](#installation)
  - [Software requirements](#software-requirements)
- [Software deployment method](#software-deployment-method)
- [Configuration](#configuration)
  - [Taxonomic assignments](#taxonomic-assignments)
  - [Preprocessing](#preprocessing)
  - [Chimera filtering](#chimera-filtering)
  - [ASV clustering tools](#asv-clustering-tools)
  - [Noise filtering](#noise-filtering)

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

or visit the [release page](https://github.com/insect-biome-atlas/happ/releases) and download the latest release.

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

#### Installation with Conda

If you prefer to use Conda, you can create a new environment with the required
software by running:

```bash
conda env create -f environment.yml
```

Then activate the environment with:
  
  ```bash
  conda activate happ
  ```

> [!TIP]
> If you are running the workflow on a high performance computing system with
> the SLURM workload manager use the configuration profile under
> `profiles/slurm`. If you're running on the Dardel HPC system at PDC
> specifically, use the configuration profile under `profiles/dardel`. See the
> README files in the respective subdirectory for more information.

## How to run

Once you have activated the software environment (either with `pixi shell` or `conda activate happ` as described above) the basic syntax to run the workflow is:

```bash
snakemake --configfile <path-to-your-configfile.yml> --profile <slurm/dardel/local> --sdm <apptainer/conda> <additional-arguments>
```

The `--configfile` argument specifies the path to a [configuration file](#configuration-file) in YAML format. The `--profile` argument specifies the [configuration profile](#configuration-file) to use. The `--sdm` flag is a short-hand for `--software-deployment-method` and specifies how Snakemake will handle [rule-specific dependencies](#software-deployment-method). See below for a description of each of these arguments.

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

To use a profile, simply add `--profile <profile-name>` to the Snakemake command line call. For example, to run the workflow on a system with SLURM you would run:

### Software deployment method

Once you have installed the main software packages needed to run the workflow
you also have to specify how Snakemake will handle rule-specific software
dependencies. This is controlled by the `--sdm` command line argument (short-hand for `--software-deployment-method`). We recommend to use [Apptainer](https://apptainer.org/) if this is available on your system. To use apptainer, add `--sdm apptainer` to the Snakemake command line call. Alternatively you can use [Conda](https://docs.conda.io/en/latest/) in which case you would use `--sdm conda`.

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
before running the ASV clustering tools. Splitting the data means that ASVs that do not share the same `split_rank` are not compared for clustering which means that you should set this parameter to a relatively high rank in the taxonomy tree. The default is `Family`.

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

#### VSEARCH assignments

To run the `vsearch` taxonomic assignment method (implemented in QIIME2) you need to provide paths to files compatible with QIIME2 under the `qiime2:` config section. 

- `ref:` This is the path to a fasta file containing reference sequences for the taxonomic assignment with the QIIME2 `feature-classifier` plugin.
- `taxfile:` This is the path to a file containing the taxonomy for the reference sequences in the fasta file.
- `ranks:` This is a list of taxonomic ranks present in the reference.

See the [QIIME2 docs](https://docs.qiime2.org/2024.10/tutorials/feature-classifier/) for how to obtain and format these files.

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