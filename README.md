[![Pixi Badge](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/prefix-dev/pixi/main/assets/badge/v0.json)](https://pixi.sh)

# ASV-clustering

## Overview

This repository contains a Snakemake workflow to filter and cluster Amplicon
Sequence Variants (ASVs). 

Currently the following clustering tools are supported:

| Software  | Reference                                                                                      | Code                                                                    |
|-----------|------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------|
| SWARM     | [Mahé et al 2014](https://peerj.com/articles/593/)                                             | [GitHub](https://github.com/torognes/swarm)                             |
| OptiClust | [Westcott & Schloss 2017](https://journals.asm.org/doi/10.1128/mSphereDirect.00073-17)         | [GitHub](https://github.com/SchlossLab/Westcott_OptiClust_mSphere_2017) |
| dbOTU3    | [Olesen et al 2017](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0176335) | [GitHub](https://github.com/swo/dbotu3)                                 |

The idea with this workflow is to make it easy to run OTU clustering with many 
different parameter settings then evaluate which settings you think works best
for your data. `vsearch` makes up the basis of the workflow by creating pairwise
alignments of the sequences (often the most time-consuming step) and several 
clustering runs can then be executed without having to re-create the alignments
for each run.

## Requirements

This workflow was developed for a Linux system. The easiest way to get set up
is to use the [pixi](https://pixi.sh/latest/) package manager to handle all 
software requirements.

## Installing

1. Clone the repository and change into the directory:

```bash
git clone git@github.com:johnne/ASV-clustering.git
cd ASV-clustering
```

2. Install pixi

If you don't have [pixi](https://pixi.sh/) installed, follow the
[instructions](https://pixi.sh/latest/#installation) to install it on your
system. Once pixi is installed you're ready to start using the workflow.

## Running the workflow

This workflow was written in
[Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) but is
executed using pixi via the `pixi.toml` file in the root of the repository. The
basic commandline syntax for running the workflow is:

```bash
pixi run <task> --configfile <path-to-your-config-file> <additional snakemake options>
```

The `tasks` configured in the `pixi.toml` file are:

- `chimera_filtering`: identifies and removes chimeras in the ASV dataset
- `clustering`: clusters the ASVs using the tools specified in the configuration file
- `numts_filtering`: identifies and removes ASV clusters that are likely to derive from Nuclear mitochondrial DNA (NUMTs)
- `all`: runs all the steps above

The config file should be a file in YAML-format specifying various parameters,
most importantly the path to your input dataset. Read more about how to
configure the workflow under [Configuration](#configuration).

The additional options are any [command line
options](https://snakemake.readthedocs.io/en/stable/executing/cli.html#all-options)
recognized by Snakemake. The most important options are already pre-configured
in the [Configuration profiles](#configuration-profiles).

### Basic usage

The workflow is preconfigured to run on a small test dataset, so you can try it
out directly by running:

```bash
pixi run all -n
```

Here `all` is the name of a task which triggers snakemake to run the entire
workflow from start to finish. The `-n` option is passed to Snakemake and means
that a 'dry-run' is performed which will give you a summary of jobs without
doing anything. Remove the `-n` flag to actually run the steps.

You will see that pixi installs the required software such as `snakemake`,
`biopython` *etc.* as well as a number of other dependencies. You will then get
a summary of the jobs that would be run.

Similarly, to get a summary of steps in the first `chimera_filtering` part of
the workflow run:

```bash
pixi run chimera_filtering -n
```

And for the `clustering` part:

```bash
pixi run clustering -n
```

And finally, for the `numts_filtering` part:

```bash
pixi run numts_filtering -n
```

The `clustering` part depends on output from the `chimera_filtering`, and
`numts_filtering` depends on output from `clustering`. Snakemake keeps track of
which steps need to be run, so you can in principle run the complete workflow in
one go with `pixi run all`. **However**, we highly recommend that you run each
step in succession, starting with `chimera_filtering`, then `clustering` and
finally `numts_filtering`. Not only does this give you a chance to inspect the
intermediate output, but also it allows the workflow to update the list of taxa
to process in subsequent steps.

### Configuration profiles

The workflow comes with three different configuration profiles. These profiles
set a number of command line options for Snakemake and define resource usage for
jobs. These can be used by adding `--profile <profile name>` to the commandline.

| Profile name | Description | Command line option |
| ------------ | ----------- | ------------------- |
| local | for running the workflow locally on your own computer | `--profile local` |
| dardel | for running on the high performance computing system [Dardel](https://www.pdc.kth.se/hpc-services/computing-systems/about-the-dardel-hpc-system-1.1053338) | `--profile dardel` |
| slurm | a general purpose profile for running on HPC systems other than Dardel | `--profile slurm` |

## Configuration

The workflow can be configured using a configuration file in YAML format.
A default config file is available at [config/config.yml](config/config.yml).
It's recommended that you copy this file and make your changes in the copy, then
pass it on the command line with `--configfile <path-to-your-configfile>`.

The configurable parameters in the config file are explained below.

### Input parameters

The `rundir` parameter is the name of a subdirectory under data/ that **must** contain:

- asv_seqs.fasta (ASV sequences in FASTA format) 
- asv_counts.tsv (tab separated file with counts of ASVs (rows) in samples (columns)) 
- and asv_taxa.tsv (tab separated file taxonomic assignments of each ASV)

The `run_name` parameter designates a name of the workflow run and can be used
to separate runs with different parameters of the clustering tools. This allows
you to try different settings of the clustering tools without having to rerun
the entire workflow from start to finish.

> [!IMPORTANT]
>
>The top-level `run_name` parameter setting is tightly linked to the `run_name` 
>config parameter under `chimera`. If you change the chimera filtering settings, 
> be sure to also update the chimera `run_name`.

The `split_rank` parameter is used to split the input ASVs by a taxonomic rank
prior to clustering. For example, setting `split_rank: "Family"` (the default)
splits the ASVs by family assignments (as given in the `asv_taxa.tsv` file).
Each split of ASVs can then be clustered in parallell which can help speed
things up especially on large datasets. The workflow will by default run on all
unique taxa names at rank `split_rank` found in the `data/{rundir}/asv_taxa.tsv`
file, and will list the taxa in the file `data/{rundir}/{split_rank}.txt`.

The `ranks` parameter lists the ranks given in the `asv_taxonomy.tsv` file and
determines the columns reported in the `cluster_consensus_taxonomy.tsv` results
file.

### Taxonomic parameters

The `evaluation_rank` parameter specifies a taxonomic rank in the
`asv_taxonomy.tsv` which will be considered as "ground truth" when calculating
precision and recall values for the clustering. By default, this is set to
"Species" which means that the precision and recall values will reflect how
often the generated clusters only contain ASVs assigned to a single species, and
how often ASVs from the same species are clustered together, respectively.

The `consensus_ranks` parameters specifies what ranks to use when attempting to
assign a conensus taxonomy to generated clusters. The default is `["Family",
"Genus", "Species"]`.

> [!NOTE]
>
>Please note that column matching between the `asv_taxa.tsv` file and the
>`split_rank`, `evaluation_rank` and `ranks` parameters is case-sensitive and
>that the `asv_taxa.tsv` file **must** have ranks in the header that match with
>your configuration settings.

The `consensus_threshold` parameter is used when assigning consensus taxonomies
to clusters and is the threshold (in %) at which the abundance weighted
taxonomic assignments for ASVs in a cluster must agree in order to assign the
taxonomy to the cluster. As an example, if a cluster contains 3 ASVs with the
following taxonomic assignments and total sum of counts across samples:

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

### Chimera filtering parameters

> [!IMPORTANT] 
>
> Running the chimera filtering on your input data may in some cases result in
>taxa with 0 sequences. The workflow cannot take this into account and will fail
>in downstream steps. We therefore recommend to **first** run the chimera
>filtering step by using `chimera_filtering` as a target to snakemake, _e.g._:
>
>```commandline
>pixi run chimera_filtering
>```
>
>Once this part of the workflow is done, on subsequent runs of the workflow only
>taxa with ASVs remaining after chimera filtering will be used as input to
>clustering.

The workflow supports optional chimera removal using the uchime algorithm 
implemented in vsearch. Chimera detection can be run either in 'batchwise' mode
using the data under `rundir` directly or in 'samplewise' mode in which the 
`asv_counts.tsv` file is used to generate sample-specific fasta files containing 
sequences with a count > 0 in each sample. Sequences identified as non-chimeric
are then passed downstream in the workflow and used as input for clustering.

Config parameters for the chimera removal part of the workflow are nested 
under `chimera:` in the config file and below are the default values:

```yaml
chimera:
  run_name: "chimera1"
  remove_chimeras: True
  method: "samplewise"
  algorithm: "uchime_denovo"
  min_samples_shared: 1
  min_frac_samples_shared: 0.5
  min_chimeric_samples: 0
  min_frac_chimeric_samples: 0
  dn: 1.4
  mindiffs: 3
  mindiv: 0.8
  minh: 0.28
```

The `run_name` parameter works in a similar manner as the top level `run_name`
and allows you to define different runs of the chimera filtering. A fasta file
with non-chimeric sequences will be produced under
`results/chimera/<rundir>/filterred/<run_name>/<method>.<algorithm>/nonchimeras.fasta`. 
This fasta file will be used as input to the rest of the workflow.

The `remove_chimeras` parameter is a boolean controlling whether to run chimera
filtering or not. Set this to False to skip chimera filtering.

The `method` parameter can be either `batchwise` (all ASVs are taken into
consideration) or `samplewise` (only ASVs occurring in the same samples are
considered).

The `algorithm` parameter can be either `uchime_denovo`, `uchime2_denovo` or
`uchime3_denovo`. The latter two require perfect matches between the ASV
sequence and a chimeric model, whereas 'uchime_denovo' does not.

The `min_samples_shared` parameter is specific to the batchwise mode and
specifies the number of samples in which chimeric ASVs have to be present with
their 'parents' (see Uchime
[docs](https://drive5.com/usearch/manual/chimeras.html)) in ordered to be marked
as chimeric.

The `min_frac_samples_shared` parameter is similar to `min_samples_shared`, but
instead of an absolute number require that sequences are present with their
parents in a fraction of the samples in which they are present.

The `min_chimeric_samples` parameter refers to the samplewise chimera mode, and
determines in how many samples an ASVs has to be marked as chimeras in to be
filtered as chimeric. Setting this value to 0 (the default) means that sequences
have to be marked as chimeric in **all** samples where they are present in order
to be filtered as chimeric.

The `min_frac_chimeric_samples` is analogous to `min_chimeric_samples` but specifies
a fraction of samples, instead of an absolute number

The `dn`, `mindiffs`, `mindiv` and `minh` parameters are specific to how the
chimeric score of sequences is calculated. Please see the Uchime
[docs](https://www.drive5.com/usearch/manual6/UCHIME_score.html) for details. 

In all algorithms ASV sequences are first aligned to other ASVs that are above a 
certain abundance threshold. This so called [abundance skew](https://www.drive5.com/usearch/manual6/abundance_skew.html)
threshold is by default set to 2.0 for the 'uchime_denovo' and 'uchime2_denovo'
algorithms and to 16.0 for 'uchime3_denovo'. However, this value can be overridden
by explicitly setting `abskew` in the config file under `chimera:`.

### Clustering tools

The workflow can be configured to run one or more ASV clustering tools. The `software` parameter takes a list of supported tool names: `swarm`, `opticlust` and `dbotu3` currently.

### Tool-specific parameters

Below are specific parameters for tools used in the workflow.

#### vsearch

Pairwise identities between ASVs are obtained by running vsearch with the
`usearch_global` setting. This is used as input by `opticlust` for clustering
ASVs.

The `threads` parameter specifies how many threads to use for the vsearch alignment.

The `id` parameter sets the minimum pairwise identity to report. The default is 0.84 meaning that sequences with a pairwise identity less than 84% is rejected.

The `iddef` parameter sets the identity definition used by vsearch. The default
setting of `1` means that identity is calculated as `(matching columns) /(alignment length)`. Refer to the vsearch manual for other settings.

The `query_cov` parameter sets the minimum aligned fraction of the query
sequence. The default setting of `0.9` means that alignments are rejected if
this fraction is less than 905.

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
[manual](https://mothur.org/wiki/cluster/).

The `aligner` parameter specifies whether to use `vsearch` (default) or the
built-in aligner in `mothur` for generating pairwise identities. Note that the
latter has not been tested for this workflow.

The `delta` parameter sets the stable value for the metric in opticlust.

The `cutoff` parameter controls at what similarity threshold clusters are
generated. 

The `initialize` parameter sets the initial randomization for opticlust. The
default is `singleton` where each sequence is randomly assigned to its own OTU.
The other accepted setting `oneotu` means that all sequences are assigned to one
otu. 

The `precision` parameter sets the floating point precision for opticlust.

The `threads` parameter specifies how many threads to allocate for the opticlust command.

#### Swarm

Swarm clusters sequences using a local linking threshold `d` which represents
the maximum number of differences between two amplicons. The input is a single
fasta sequence where the total count of each sequence is suffixed to the
sequence ID.

For a full description of Swarm parameters, see the
[manual](https://github.com/torognes/swarm/blob/master/man/swarm_manual.pdf).

The `differences` parameter specifies the maximum number of differences
(mismatches) allowed between two ASVs, meaning that two ASVs will be grouped if
they have integer (or less) differences.

The `no-otu-breaking` parameter deactivates Swarms built-in cluster refinement
when using `d=1`.

The `fastidious` parameter specifies whether to perform a second clustering pass
(when set to True) to reduce the number of small clusters when using `d=1`.

The `boundary` parameter defines the minimum abundance of what should be
considered a _large_ cluster when using `d=1`.

The `threads` parameter specifies how many threads to allocate to the swarm
command.

The parameters `match-reward`, `mismatch-penalty`, `gap-opening-penalty` and
`gap-extension-penalty` allows control of Swarms advanced Pairwise alignment
options.

#### dbotu3

dbotu3 uses both sequences and their distribution across samples to cluster ASVs
into OTUs. It takes as input ASV sequences in fasta format and a table of
counts.

For a full description of dbotu3 parameters, see the
[documentation](https://dbotu3.readthedocs.io/en/latest/).

The `dist` parameter sets the maximum allowed genetic dissimilarity between
sequences.

The `abund` parameter sets the minimum fold difference for comparing two OTUs.

The `pval` parameter sets the minimum p-value for merging ASVs.


### NUMTs filtering

The generated clusters are used as input for removal of Nuclear mitochondrial
DNA segments (NUMTs). In this part of the workflow, the representative sequences
of clusters are compared based on abundance and sequence similarity to other
clusters in the same taxonomic order to identify and remove sequences that are
likely to represent NUMTs. As an additional step, clusters that with no
taxonomic assignments at the order level (configurable with the
`filter_unclassified_rank` parameter, see below) are also removed. The original
output (prior to NUMTs removal) is kept intact however, allowing you to compare
the output at the different steps. 

The following configuration parameters determine the behaviour of the NUMTs filtering:

The `n_closest` parameter sets the number of other most closely related
sequences to compare each cluster to. A higher `n_closest` value increase the
chance of identifying NUMTs, but lead to longer runtimes.

The `non_numt_ASVs` **optional** parameter may point to a file listing so called
'trusted' ASV ids (one per line) in the dataset that are known beforehand to not
represent NUMTs. This is only used for evaluation purposes and does not
influence the actual filtering steps.

The `spikein_file` **optional** parameter may point to a tab-separated file with
taxonomic information of biological spike in taxa that have been used in the
project. If given, this file is used to identify clusters corresponding to the
spike in taxa. Those clusters are then used to calibrate the count data during
the NUMTs filtering steps. This file must have either a column named `Species`
representing species names of the spike ins, and/or a column named `BOLD_bin`
with BOLD bin ids of the spike in taxa. For the BOLD_bin column, multiple
BOLD_bin ids may be given on a single line and should then be quoted and
separated by `;`. An example of such a file is shown below:

| Species | BOLD_bin |
| ------- | -------- |
| Gryllodes sigillatus | BOLD:ABW5620 |
| Gryllus bimaculatus  | BOLD:ABX6319 |
| Shelfordella lateralis | BOLD:AAG9959 |
| Drosophila serrata | BOLD:AAU1484 |
| Drosophila bicornuta | BOLD:AAV6760 |
| Drosophila jambulina | "BOLD:AAV6734; BOLD:AAL1587" |

The `spikein_method` parameter determines whether the `Species` or `BOLD_bin`
column in the `spikein_file` should be used to identify clusters.

The `large_orders` **optional** parameter is a list of taxonomic orders known to represent a large fraction of the dataset, and to consequently contain a large set of ASV clusters. During NUMTs filtering a distance matrix is generated for all ASV clusters within each order. This step becomes very resource demanding with a lot of input sequences. By specifying such large orders here the workflow will allocate more memory and runtime when analysing these orders. 

The `filter_unclassified_rank` parameter specifies a taxonomic rank at which to filter out unassigned clusters. The default value `Order` means that clusters with no taxonomic assignment at Order level are removed as part of the NUMTs filtering steps. Note that when this parameter is set to `order` (case insensitive), NUMTs filtering will not be performed for unassigned orders (*i.e.* orders starting with `unclassified`).

### Workflow output

The final output of the workflow is placed in subdirectories for each tool used. The path to the results is determined by your configuration settings.

Below are some examples of output obtained using the default [config/config.yml](config/config.yml).

#### Chimera filtering output

The default config file has the following settings for input and chimera filtering:

```yaml
rundir: "test"
run_name: "run1"
split_rank: "Family"
chimera:
  run_name: "chimera1"
  method: "samplewise"
  algorithm: "uchime_denovo"
```

Running the following on the command line:

```bash
pixi run chimera_filtering --configfile config/config.yml
```

would result in

```
results/chimera/test/filtered/chimera1/samplewise.uchime_denovo/
├── nonchimeras.fasta # ASV sequences identified as non-chimeric (used in clustering)
├── Family.txt # List of families remaining after chimera filtering (used in clustering)
├── orders.txt # List of orders remaining after chimera filtering (used in numts filtering)
└── chimeras.fasta # ASV sequences identified as chimeric
```

So here the output follows the generalized pattern:

```
results/chimera/{rundir}/filtered/{chimera_run}/{method}.{algorithm}/
```

where 

- `rundir` corresponds to the config parameter with the same name
- `chimera_run` is the `run_name` config parameter nested under the `chimera` config entry
- `method` is the method ('samplewise' or 'batchwise') chimera config parameter
- `algorithm` is the vsearch algorith used (`uchime_denovo`, `uchime2_denovo` or `uchime3_denovo`)


#### Clustering output

The `nonchimeras.fasta` file generated in the chimera filtering is passed on to
the clustering tools defined in the config file. The default config file has

```yaml
software:
  - "swarm"
```

which means that only swarm will be run on the chimera-filtered sequences.

Running:

```bash
pixi run clustering --configfile config/config.yml
```

would result in the following:

```
results/swarm/test/chimera1/samplewise.uchime_denovo/Family/runs/run1/
├── cluster_taxonomy.tsv
├── cluster_reps.fasta
├── cluster_counts.tsv
├── cluster_consensus_taxonomy.tsv
├── precision_recall.txt
└── precision_recall.order.txt
```

where

**cluster_taxonomy.tsv** is a tab-separated file of all non-chimeric ASVs (rows) with the following columns:

| column name | description |
| ----------- | ----------- |
| ASV         | ASV id      |
| cluster     | name of the cluster that the ASV belongs to |
| median      | the median relative abundance of the ASV in the dataset |
| mean        | the mean relative abundance of the ASV in the dataset |
| sum         | the summed relative abundance of the ASV in the dataset |
| rank1...rankN | columns matching the ranks supplied in the input `asv_taxa.tsv` file |
| representative | 1 if the ASV was selected as representative of the cluster, otherwise 0

The `cluster` name is prefixed with the corresponding taxonomic label assigned
to the `split_rank` rank used. So with the default config which has `split_rank:
Family` the clusters would be prefixed with the assigned Family name of the ASV.

**cluster_reps.fasta** is a fasta file of sequences for all representative ASVs,
the headers contain both the ASV id and the corresponding cluster name.

**cluster_counts.tsv** is a tab-separated file with the abundance of clusters
(rows) in samples (columns) obtained by summing the counts of all ASVs contained
in each cluster.

**cluster_consensus_taxonomy.tsv** is a tab-separated file where the taxonomy of
clusters have been resolved using a consensus approach (described above) that
takes into account the taxonomic assignments of all ASVs within a cluster.

**precision_recall.txt** is a text file containing statistics on the clustering.
This file shows how well the clustering output corresponds to the level of
taxonomy specified by the `evaluation_rank` in the config file. In the default
config file, `evaluation_rank` is set to 'Species' which means that the
statistics will show how the generated clusters correspond to species level
assignments. For the purpose of this evaluation, only ASVs classified with
unambiguous assignments at `evaluation_rank` are used. The statistics shown are:

| statistic | description |
| --------- | ----------- |
| Total clusters | Number of total clusters evaluated (after removing unclassified and ambiguous assignments) |
| Total number of `<evaluation_rank>` | Number of unique taxonomic labels at the level of `evaluation_rank` |
| Total positives | Total number of pairwise ASV comparisons analysed | 
| True positives | Number of times two compared ASVs from the same `evaluation_rank` were found in *the same* cluster |
| False positives | Number of times two compared ASVs from the same `evaluation_rank` were found in *different* clusters |
| False negatives | Number of times two compared ASVs from *different* `evaluation_rank` were found in the same cluster |
| precision | A measure of how often clusters contain only one unique `evaluation_rank` |
| recall | A measure of how often ASVs from the same `evaluation_rank` are placed in the same cluster |
| homogeneity/completeness | These values are calculated using the [homogeneity_score](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.homogeneity_score.html) and [completeness_score](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.completeness_score.html) in `scikit-learn`. They mean essentially the same thing as precision/recall but are not as fine-grained on the clustering output. |

So if we are evaluating the clusters against species level assignments, then
`precision` indicates to what degree the generated clusters contain only one
unique species, while `recall` indicates how well the workflow placed ASVs from
the same species into the same cluster. It is possible to maximize the precision
by under-clustering the data, *i.e.* generating a large number of 'pure'
clusters which means that most ASVs from the same species will be in different
clusters. And similarly, recall can be maximized by over-clustering, *i.e.*
generating large 'impure' clusters with ASVs from several unique species. In our
experience, the supported tools are more prone to under-clustering with regards
to species (meaning they give high precision values). Typically a precision
value >0.9 and a recall value >0.8 should be obtained, but this may vary
depending on your dataset (and taxonomic assignments).

**precision_recall.order.txt** is a text file with the same type of information as `precision_recall.txt` but with results shown for each taxonomic orders

#### Statistics summary output

The clustering statistics described above are also summarized for each clustering tool used under `results/stats/{rundir}/{chimera_run}/{method}.{algorithm}/{rank}/runs/{run_name}.tsv` which allows for quick comparison of the results from several tools.

#### NUMTs filtering output

The NUMTs filtering part of the workflow generates three files in the same output directory as the clustering, steps:

- **non_numts.tsv** is a tab-separated with the same columns as the **cluster_taxonomy.tsv** file described above under [Clustering output](#clustering-output), but only containing ASVs in non-numts cluster.

- **non_numts_clusters.fasta** contains the representative ASV sequences of clusters remaining after NUMTS filtering.

- **precision_recall.non_numts.txt** gives the same type of statistics on the clustering, but only taking into account the clusters remaining after NUMTs filtering.

In addition, a subfolder called `numts_filtering` contains the intermediate output as well as evaluation files generated during the filtering steps. This folder contains output for each taxonomic order:

```
numts_filtering/
├── aa/ # protein translated sequences
├── abundance_filter/ # output from filtering based on cluster abundance
├── alignments/ # output from PAL2NAL
├── cluster_analysis/ # intermediate summary files
├── combined_filter/ # output from combined filtering using both abundance and sequence similarity
├── evaluation/ # evaluation statistics output
├── mafft/ # output from mafft alignment
└── trimmed/ # trimmed sequences
```