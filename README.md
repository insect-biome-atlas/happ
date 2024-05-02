# ASV-clustering

## Overview

This repository contains a Snakemake workflow for clustering Amplicon Sequence Variants (ASVs). Currently the following tools are supported:

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

- Linux system
- [pixi](https://prefix.dev/) software manager

## Installing

Clone the repository and change into the directory:

```bash
git clone git@github.com:johnne/ASV-clustering.git
cd ASV-clustering
```

If you don't have [pixi](https://prefix.dev/) installed, install it by running:

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

## Quickstart

To try the workflow on a small test dataset, simply run:

```bash
pixi run local
```

## Detailed instructions

The workflow is configured using the file [config/config.yml](config/config.yml). To modify the workflow behaviour simply change settings directly in this file.

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

### Chimera detection

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

The `run_name` parameter works in a similar manner as the base level `run_name`
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
specifies the number of  samples in which chimeric ASVs have to be present with
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

### Aligning

Pairwise identities between ASVs are obtained by running vsearch with the
`usearch_global` setting. The `vsearch` config entry determines the behaviour of this step.

The `threads` parameter specifies how many threads to use for the vsearch alignment.

The `id` parameter sets the minimum pairwise identity to report. The default is 0.84 meaning that sequences with a pairwise identity less than 84% is rejected.

The `iddef` parameter sets the identity definition used by vsearch. The default
setting of `1` means that identity is calculated as `(matching columns) /(alignment length)`. Refer to the vsearch manual for other settings.

The `query_cov` parameter sets the minimum aligned fraction of the query
sequence. The default setting of `0.9` means that alignments are rejected if
this fraction is less than 905.

### ASV clustering

The workflow can be configured to run one or more ASV clustering tools. The `software` parameter takes a list of supported tool names: `swarm`, `opticlust` and `dbotu3` currently.

- opticlust
- swarm
- dbotu3

Below are descriptions of the configurable tool-specific parameters.

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

The `delta` parameter sets the stable value for the metric in opticlust Default
delta=0.0001. 

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


### Workflow output

The final output of the workflow is placed in subdirectories for each tool used. The path to the results is determined by your configuration settings and follows the generalized pattern:

`results/{tool}/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/`

where 

- `tool` can be `swarm`, `opticlust` or `dbotu3`
- `rundir` corresponds to the config parameter with the same name
- `chimera_run` is the `run_name` config parameter nested under the `chimera` config entry
- `chimdir` is a combination of the `method` ('samplewise' or 'batchwise') and `algorithm` config parameter
- `rank` is the config setting for `split_rank`
- `run_name` is the top-level config parameter with the same name

As an example, with the following config settings:

```yaml
rundir: "project1"
run_name: "run1"
split_rank: "Family"
software: ["swarm", "opticlust"]
chimera:
  run_name: "sw.strict"
  method: "samplewise"
  algorithm: "uchime_denovo"
```

there would be results created under:

`results/swarm/project1/sw.strict/samplewise.uchime_denovo/Family/runs/run1/`

and

`results/opticlust/project1/sw.strict/samplewise.uchime_denovo/Family/runs/run1/`

The following files are created by the workflow:

**cluster_taxonomy.tsv**

This is a tab-separated file containing all (filtered) ASVs (rows) along with the following columns:

- cluster: Name of the ASV cluster that the ASV belongs to
- median: The median of normalized counts of the ASV across all samples
- mean: The mean of normalized counts of the ASV
- sum: The sum of normalized counts of the ASV
- <taxonomic ranks>: The taxonomic assignments taken from the `asv_taxa.tsv` input file
- representative: A column with value = 1 if the ASV is a representative of the cluster

**cluster_consensus_taxonomy.tsv**

Assigned consensus taxonomy of clusters.

**cluster_counts.tsv**

Counts of clusters (rows) across samples (columns). The counts of clusters are
generated by summing counts for all ASVs in the cluster.

**cluster_reps.fasta**

Sequences of cluster representatives in FASTA format.

#### Stats output

The workflow will evaluate the clustering using the taxonomic assignments in
`data/<rundir>/asv_taxa.tsv` at the `evaluation_rank` set in the configuration
file as the 'ground truth' and calculate various statistics such as
precision/recall, completeness/homogeneity of the clusters. 

- Precision is calculated as: `precision = TP/(TP + FP)`
- Recall is calculated as: `recall = TP/(TP + FN)`

where

- True positives (TP) is calculated by counting the number of times two ASVs belonging 
  to the same ASV cluster also share the same taxa at `evaluation_rank`.
- False positives (FP) is calculated as: `FP = totalPositives - TP`.
- False negatives (FN) is calculated by counting the number of times ASVs that share
  the same taxa at `evaluation_rank` are placed into different ASV clusters.

The evaluation script used in this workflow begins by calculating the total number
of positives from a given dataframe. This is formulated as: `totalPositives = N * (N - 1) / 2`
where `N` is the total number of ASVs.

Statistics are placed under `results/stats/{rundir}/{chimera_run}/{chimdir}/{rank}/runs/{run_name}/` so if using the example configuration above there would be: `results/stats/project1/sw.strict/samplewise.uchime_denovo/Family/runs/run1.tsv`

Where the file `run1.tsv` summarizes the evaluation of all tools used. 

Note that when calculating these statistics, sequences unassigned at
`evaluation_rank` (Species by default) are ignored, so the total number of ASVs,
clusters and species is likely higher.

### 7. Evaluation
True and false positives, as well as true and false negatives can be evaluated
for each clustering results using the taxonomic assignments of ASVs. The `evaluation_rank`
config parameter (default = "Species") determines what to evaluate the clustering
against (whatever assignment is given for this rank that is taken as the ground
truth for an ASV). 


