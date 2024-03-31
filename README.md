# ASV-clustering

## Overview

This repository contains a Snakemake workflow for clustering Amplicon Sequence Variants (ASVs). Currently the following tools are supported:

| Software  | Reference                                                                                      | Code                                                                    |
|-----------|------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------|
| SWARM     | [Mahé et al 2014](https://peerj.com/articles/593/)                                             | [GitHub](https://github.com/torognes/swarm)                             |
| OptiClust | [Westcott & Schloss 2017](https://journals.asm.org/doi/10.1128/mSphereDirect.00073-17)         | [GitHub](https://github.com/SchlossLab/Westcott_OptiClust_mSphere_2017) |
| dbOTU3    | [Olesen et al 2017](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0176335) | [GitHub](https://github.com/swo/dbotu3)                                 |

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

The `run_name` parameter designates a name of the workflow run and can be used to separate runs with different parameters of the clustering tools. This allows you to try different settings of the clustering tools without having to rerun the entire workflow from start to finish.

The `split_rank` parameter is used to split the input ASVs by a taxonomic rank prior to clustering. For example, setting `split_rank: "Family"` (the default) splits the ASVs by family assignments (as given in the `asv_taxa.tsv` file). Each split of ASVs can then be clustered in parallell which can help speed things up especially on large datasets. The final output of the workflow will still be combined ASV cluster files.

The `ranks` parameter lists the ranks given in the `asv_taxonomy.tsv` file and determines the columns reported in the `cluster_consensus_taxonomy.tsv` results file.

The `consensus_ranks` parameters specifies what ranks to use when attempting to assign a conensus taxonomy to generated clusters. The default is `["Family", "Genus", "Species"]`.

> [!NOTE]
>
>Please note that column matching between the `asv_taxa.tsv` file and the `split_rank`, `evaluation_rank` and `ranks` parameters is case-sensitive and that the `asv_taxa.tsv` file **must** have ranks in the header that match with your configuration settings.

The `consensus_threshold` parameter is used when assigning consensus taxonomies to clusters and 
is the threshold (in %) at which the abundance weighted taxonomic assignments for ASVs in a 
cluster must agree in order to assign the taxonomy to the cluster.
As an example, if a cluster contains 3 ASVs with the following taxonomic assignments and total 
sum of counts across samples:

| ASV | Family | Genus | Species | ASV_sum |
|-----|--------|-------|---------|---------|
| ASV1 | Tenthredinidae | Pristiphora | Pristiphora mollis | 20 |
| ASV2 | Tenthredinidae	| Pristiphora	| Pristiphora cincta | 60 |
| ASV3 | Tenthredinidae	| Pristiphora	| Pristiphora leucopodia | 20 |

then at Species the abundance weighted taxonomic assignment is 60% _Pristiphora cincta_ and 20% 
_Pristiphora mollis_ and _Pristiphora leucopodia_ each. At a 80% consensus threshold we cannot 
assign a taxonomy at Species level to the cluster, so the algorithm will move up to Genus level 
where we have 100% agreement for Pristiphora. The final consensus taxonomy for the cluster will 
then be:

| Cluster | Family | Genus | Species |
|---------|--------|-------|---------|
| Cluster1 | Tenthredinidae | Pristiphora | unresolved.Pristiphora |

At `consensus_threshold: 60` we would have been able to assign taxonomy at the Species level and 
the cluster taxonomy would have been:

| Cluster | Family | Genus | Species |
|---------|--------|-------|---------|
| Cluster1 | Tenthredinidae | Pristiphora | Pristiphora cincta |

## Workflow overview
The idea with this workflow is to make it easy to run OTU clustering with many 
different parameter settings then evaluate which settings you think works best
for your data. `vsearch` makes up the basis of the workflow by creating pairwise
alignments of the sequences (often the most time-consuming step) and several 
clustering runs can then be executed without having to re-create the alignments
for each run.

### 0. Chimera detection (optional)

> [!IMPORTANT]
>Please note that running the chimera filtering on your input data may in some
>cases result in taxa with 0 sequences. The workflow cannot take this into account
>and will fail in downstream steps. We therefore recommend to **first** run the chimera
>filtering step by using `chimera_filtering` as a target to snakemake, _e.g._:
>
>```commandline
>snakemake --profile local --configfile <yourconfig.yaml> -c 4 chimera_filtering
>```
>
>Once this part of the workflow is done, on subsequent runs of the workflow only taxa
>with ASVs remaining after chimera filtering will be used as input to clustering.

The workflow supports optional chimera removal using the uchime algorithm 
implemented in vsearch. Chimera detection can be run either in 'batchwise' mode
using the data under `rundir` directly or in 'samplewise' mode in which the 
`asv_counts.tsv` file is used to generate sample-specific fasta files containing 
sequences with a count > 0 in each sample. Sequences identified as non-chimeric
are then passed downstream in the workflow and used as input for clustering.

#### Parameters for chimera detection and filtering
Config parameters for the chimera removal part of the workflow are nested 
under `chimera:` in the config file and below are the default values.

```yaml
chimera:
  run_name: "chimera1"
  remove_chimeras: True
  method: "samplewise"
  algorithm: "uchime_denovo"
  min_samples_shared: 1
  min_frac_samples_shared: 0.5
  min_chimeric_samples: 0
  dn: 1.4
  mindiffs: 3
  mindiv: 0.8
  minh: 0.28
```

- `run_name`: this parameter allows you to define different runs of the 
chimera filtering. A fasta file with non-chimeric sequences will be produced 
under `results/chimera/<rundir>/filterred/<run_name>/<method>.<algorithm>/nonchimeras.fasta`.
- `remove_chimeras`: set this to False to skip chimera filtering
- `method`: run chimera detection in `batchwise` or `samplewise` mode.
- `algorithm`: select between `uchime_denovo`, `uchime2_denovo` or 
  `uchime3_denovo`. The latter two require perfect matches between the ASV
  sequence and a chimeric model, whereas 'uchime_denovo' does not.
- `min_samples_shared`: in batchwise mode, chimeric ASVs have to be present 
  with their 'parents' (see Uchime [docs](https://drive5.com/usearch/manual/chimeras.html)) 
  in `min_samples_shared` number of samples in order to be filtered as chimeric.
- `min_frac_samples_shared`: similar to `min_samples_shared`, but instead of 
  an absolute number require that sequences are present with their parents 
  in a fraction of the samples
- `min_chimeric_samples`: in samplewise mode, chimeric ASVs have to be 
  marked as chimeras in `min_chimeric_samples` number of samples in order to 
  be filtered as chimeric. Setting this value to 0 (the default) means that 
  sequences have to be marked as chimeric in **all** samples where they are 
  present in order to be filtered as chimeric.

The `dn`, `mindiffs`, `mindiv` and `minh` parameters are specific to how the 
chimeric score of sequences is calculated. Please see the Uchime [docs](https://www.drive5.com/usearch/manual6/UCHIME_score.html) 
for details. 

In all algorithms ASV sequences are first aligned to other ASVs that are above a 
certain abundance threshold. This so called [abundance skew](https://www.drive5.com/usearch/manual6/abundance_skew.html)
threshold is by default set to 2.0 for the 'uchime_denovo' and 'uchime2_denovo'
algorithms and to 16.0 for 'uchime3_denovo'. However, this value can be overridden
by explicitly setting `abskew` in the config file under `chimera:`, _e.g._:

The `run_name` parameter nested under `chimera:` in the config file
controls the output structure of the chimeric filtering part of the workflow and
influences input to downstream rules. Changing this value allows you to 
tune the chimeric detection without rerunning unnecessary jobs. 

For example:

Say you have input data under `data/project1/` with counts of ASVs in two samples
`sample1` and `sample2`. Then the `rundir` config parameter should be set to 
`project1`. 

Now say you want to run chimera detection using the `batchwise` method and the 
`uchime_denovo` algorithm. Then you should set:

```yaml
chimera:
  remove_chimeras: True
  method: "batchwise"
  algorithm: "uchime_denovo"
```
in your config file. In addition, setting the `run_name` under `chimera` in your
config file gives the chimera detection with this set of parameters a name:

```yaml
chimera:
  remove_chimeras: True
  method: "batchwise"
  algorithm: "uchime_denovo"
  run_name: "chimera1"
```

Running 

```commandline
snakemake --profile local --configfile <yourconfig.yaml> -c 1 chimera_filtering 
```

will generate:

```
results/chimera/project1
├── batchwise.uchime_denovo
│   └── uchimeout.txt
└── filtered
    └── chimera1
        └── batchwise.uchime_denovo
            ├── chimeras.fasta
            ├── nonchimeras.fasta
            └── uchimeout.tsv
```

This means that the main output from chimera detection (with `vsearch`) using the `minh`,
`mindiffs` and `mindiv` defined in the config file are placed under 
`results/chimera/project1/batchwise.uchime_denovo/`. The `uchimeout.txt` file in
this folder is then used to filter chimeras with the `min_samples_shared` and 
`min_frac_samples_shared` values set in the same config file, resulting in the 
`chimeras.fasta` and `nonchimeras.fasta` files under 
`results/chimera/project1/filtered/chimera1/batchwise.uchime_denovo/` folder.

Changing the `method` parameter to 'samplewise' and rerunning the workflow will
trigger a new run of chimera detection, this time on sample-specific fasta 
files, and result in:

```
results/chimera/project1
├── batchwise.uchime_denovo
│   └── uchimeout.txt
├── filtered
│   └── chimera1
│       ├── batchwise.uchime_denovo
│       │   ├── chimeras.fasta
│       │   ├── nonchimeras.fasta
│       │   └── uchimeout.tsv
│       └── samplewise.uchime_denovo
│           ├── chimeras.fasta
│           └── nonchimeras.fasta
└── samplewise.uchime_denovo
    └── samples
        ├── sample1
        │   └── uchimeout.txt.gz
        └── sample2
            └── uchimeout.txt.gz
```

Say that you now want to try chimera detection either with different values
of `minh`, `mindiffs` or `mindiv` or requiring that a sequence is only required
to be marked as chimeric in a single sample (controlled via the 
`min_chimeric_samples` parameter) you can modify these parameters and change the
`run_name` under the chimera config section to _e.g._ `chimera2`. Now a rerun of 
the workflow will only trigger the filtering step to be run again with the 
updated parameters. The actual vsearch chimera detection is not run again.
The resulting output will be:

```
results/chimera/project1/
├── batchwise.uchime_denovo
│   └── uchimeout.txt
├── filtered
│   ├── chimera1
│   │   ├── batchwise.uchime_denovo
│   │   │   ├── chimeras.fasta
│   │   │   ├── nonchimeras.fasta
│   │   │   └── uchimeout.tsv
│   │   └── samplewise.uchime_denovo
│   │       ├── chimeras.fasta
│   │       └── nonchimeras.fasta
│   └── chimera2
│       └── samplewise.uchime_denovo
│           ├── chimeras.fasta
│           └── nonchimeras.fasta
└── samplewise.uchime_denovo
    └── samples
        ├── sample1
        │   └── uchimeout.txt.gz
        └── sample2
            └── uchimeout.txt.gz
```

**Note that the setting for the `run_name` under the chimera parameters will
carry downstream in the workflow. Disabling chimera detection will 

### 1. Splitting the input
The workflow is set up to split the input data (either raw or chimera filtered)
by taxonomic rank family and perform clustering of sequences in each family in 
parallell. However, the data can be arbitrarily split on any taxonomic rank, 
configurable via the `split_rank` config parameter. The workflow will by default 
run on all unique taxa names at rank `split_rank` found in the 
`data/{rundir}/asv_taxa.tsv` file, and will list the taxa in the file 
`data/{rundir}/{split_rank}.txt`. You can create that latter file yourself if 
you want, which will override the default behaviour and cause 
the workflow to only run on taxa listed. For example, if you list (one per row)
only your families of interest in `data/{rundir}/Family.txt` then the workflow
will only run with sequences belonging to those families.
Note that if you do not specify a file with ranks, _e.g._ `Family.txt` then the
workflow will read the `asv_taxa.tsv` file, identify taxa at `{split_rank}` but 
only output those that have at least one ASV with a total count > 0 in the counts
file **and** that have a sequence in the `asv_seqs.fasta` file. 

> [!IMPORTANT]
>Please note that running the chimera filtering on your input data may in some
>cases result in taxa with 0 sequences. The workflow cannot take this into account
>and will fail in downstream steps. We therefore recommend to **first** run the chimera
>filtering step by using `chimera_filtering` as a target to snakemake, _e.g._:
>
>```commandline
>snakemake --profile local --configfile <yourconfig.yaml> -c 4 chimera_filtering
>```
>
>Once this part of the workflow is done, on subsequent runs of the workflow only taxa
>with ASVs remaining after chimera filtering will be used as input to clustering.

### 2. Filtering and formatting
Initially, the counts for each sequence is summed across samples and any potential 
sequence with a sum of 0 is removed. The filtering step also ensures that the 
taxonomic and sequence data lines up with no missing sequences in either. If you
have set chimera detection to True (see above) then this part of the workflow
will use the nonchimeric sequences identified.

### 3. Aligning
Filtered sequences are then aligned per `split_rank` (Family by default) using
`vsearch` like so:

```bash
vsearch --usearch_global <input> --db <input> --self --userout <output-distance> \
  -userfields query+target+id --maxaccepts 0 --maxrejects 0 --id {params.id} \
  --iddef {params.iddef}  --query_cov {params.query_cov} --threads {threads}
```

where `params.id=0.84`, `params.iddef=1`, `params.query_cov=0.9` are the default
settings (configurable via the config file):

```yaml
vsearch:
  threads: 10
  id: 0.84
  iddef: "1"
  query_cov: 0.9
```

### 4. OTU clustering
When sequences have been aligned and pairwise distances are available the workflow
continues with clustering of ASVs into OTUs using the tools listed in `config["software"]`.
By default the tools are:

- opticlust
- swarm
- dbotu3

#### 4.1 opticlust
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

#### 4.2 swarm
Swarm clusters sequences using a local linking threshold `d` (set to 1 by default)
which represents the maximum number of differences between two amplicons. The input
is a single fasta sequence where the total count of each sequence is suffixed
to the sequence ID.

#### 4.3 dbotu3
dbotu3 uses both sequences and their distribution across samples to cluster ASVs
into OTUs.

### 5. Generating OTU tables
While the output from the tools used here differ, the workflow will generate a 
standardized table called `asv_clusters.tsv` with ASVs and their corresponding 
cluster for each combination of tool and taxa. This table has the following format:

| ASV   | cluster  |
|-------|----------|
| asv1  | cluster1 |
| asv2  | cluster1 |
| asv3  | cluster2 |



### 6. Workflow output
Below is an example of the results directory after a typical run of the 
workflow. Config settings in this example are shown below:

Example config:
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

#### Vsearch output
Pairwise distance output from vsearch forms the base of the clustering. This is 
placed under `results/vsearch/<rundir>/<chimera_run_name>/<chimera_method>.<chimera_algorithm>/`.

```bash
results/vsearch/project1/sw.strict
`-- samplewise.uchime_denovo
    `-- Family                      
        `-- taxa                             
            |-- Acanthosomatidae
            |   |-- asv_seqs.dist.gz
            |   `-- asv_seqs.dist.reformat.gz
            ...
```

#### Swarm output
Output files specific to swarm are placed under 
`results/swarm/<rundir>/<chimera_run_name>/<chimera_method>.<chimera_algorithm>/`.

```bash
results/swarm/project1/sw.strict
`-- samplewise.uchime_denovo                 
    `-- Family                         
        |-- runs                     
        |    -- run1           
        |       |-- asv_reps.counts.tsv 
        |       |-- asv_reps.fasta         
        |       |-- asv_reps.taxonomy.tsv
        |       |-- precision_recall.order.txt
        |       |-- precision_recall.txt     
         -- taxa                       
            |-- Acanthosomatidae 
            |   |-- derep.txt
            |   |-- reformat.fasta.gz
            |   `-- run1             
            |       |-- asv_clusters.tsv
            |       |-- asv_reps.counts.tsv
            |       |-- asv_reps.fasta     
            |       |-- asv_reps.taxonomy.tsv
            |       |-- swarm_table.tsv      
            |       `-- swarm.txt
            ...
```

#### Opticlust output

```bash
results/opticlust/project1/sw.strict 
`-- samplewise.uchime_denovo                   
    `-- Family                          
        |-- runs                             
        |   `-- run1                          
        |       |-- asv_reps.counts.tsv                                                              
        |       |-- asv_reps.fasta             
        |       |-- asv_reps.taxonomy.tsv                                                            
        |       |-- precision_recall.BOLD_bin.order.txt
        |       |-- precision_recall.BOLD_bin.txt                                                    
        |       |-- precision_recall.order.txt 
        |       `-- precision_recall.txt
        `-- taxa                             
            |-- Acanthosomatidae        
            |   `-- run1                     
            |       |-- asv_clusters.tsv
            |       |-- asv_reps.counts.tsv
            |       |-- asv_reps.fasta
            |       |-- asv_reps.taxonomy.tsv
            |       |-- asv_seqs.opti_mcc.list
            |       |-- asv_seqs.opti_mcc.sensspec 
            |       `-- asv_seqs.opti_mcc.steps
            ... 
```

#### Stats output
The workflow will evaluate the clustering using the taxonomic assignments in 
`data/<rundir>/asv_taxa.tsv` as the 'ground truth' and calculate various 
statistics such as precision/recall, completeness/homogeneity of the 
clusters. See below under **7. Evaluation** for details. By default, the 
clustering results are evaluated against the 'Species' level in the 
taxonomic data, but this can be changed with the config parameter 
`evaluation_rank`.

```bash
results/stats/project1/sw.strict
`-- samplewise.uchime_denovo
    `-- Family
        `-- runs
            |-- run1.order.tsv
            `-- run1.tsv
```

The file `run1.tsv` summarizes the evaluation of all tools used. An example 
output is shown below:

```bash
          clusters species precision recall  homogeneity  completeness ASVs
swarm     14520    10881   0.9699    0.9142  0.9952       0.9746       224562
opticlust 19835    10881   0.9997    0.7333  0.9997       0.9398       224562
```

This shows that swarm and opticlust generated 14520 and 19835 clusters 
respectively. The total number of ASVs going into the clustering was 224562 
and there were a total of 10881 unique species. Swarm obtained precision and 
recall of values of 0.9699 and 0.9142, respectively. Opticlust had slightly 
higher precision but lower recall. This indicates that Opticlust 
slightly overclustered the input data, which is also seen in the larger number 
of clusters produced.

Note that when calculating these numbers sequences not assigned to 
`evaluation_rank` (Species by default) are ignored, so the total number of 
ASVs, clusters and species is likely higher.

The file `run1.order.tsv` contains the same information, but with values 
calculated for each taxonomic Order.

### 7. Evaluation
True and false positives, as well as true and false negatives can be evaluated
for each clustering results using the taxonomic assignments of ASVs. The `evaluation_rank`
config parameter (default = "Species") determines what to evaluate the clustering
against (whatever assignment is given for this rank that is taken as the ground
truth for an ASV). 

- Precision is calculated as: `precision = TP/(TP + FP)`
- Recall is calculated as: `recall = TP/(TP + FN)`

where

- True positives (TP) is calculated by counting the number of times two ASVs belonging 
  to the same OTU cluster also share the same taxa at `evaluation_rank`
- False positives (FP) is calculated as: `FP = totalPositives - TP`
- False negatives (FN) is calculated by counting the number of times ASVs that share
  the same taxa at `evaluation_rank` are placed into different OTU clusters

The evaluation script used in this workflow begins by calculating the total number
of positives from a given dataframe. This is formulated as: `totalPositives = N * (N - 1) / 2`
where `N` is the total number of ASVs.


## Benchmarking

### Chimera benchmark

We benchmarked the performance of the chimera detection on an initial version of the ASV dataset (636297 ASVs in 5231 samples) using opticlust (implemented in mothur v1.44.11).

The following chimera filtering runs were tested:

| Name | Mode | Settings | Config file |
|------|------|----------|-------------|
| batchwise.strict | batchwise | ASVs had to share samples with their parents in 50% of samples in which they were present | `config/AllSamples_2209.run1.yml` |
| batchwise.lenient | batchwise | ASVs had to share samples with their parents in at least 1 sample | `config/AllSamples_2209.run2.yml` |
| samplewise.strict | samplewise | ASVs had to be identified as chimeric in all samples where they were present in order to be removed | `config/AllSamples_2209.run3.yml` |
| samplewise.lenient | samplewise | ASVs had to be identified as chimeric in at least 1 sample in order to be removed | `config/AllSamples_2209.run4.yml` |

These configurations were run as:

```bash
snakemake --profile slurm --configfile <config-file>
```

In order to evaluate the chimera filtering we matched our ASVs to a set of [trusted CO1 reference sequences](http://dx.doi.org/10.5883/DS-FINPRO) and considered ASVs with perfect matches to these sequences as non-chimeric in order to estimate 'false positives', *i.e.* true biological sequences identified as chimeras. We also calculated:

- total number of ASVs removed as chimeric
- number of ASVs removed that were only found in one sample
- sum of reads for removed ASVs
- sum of reads for removed 'trusted' ASVs
- ratio of generated clusters/species assignments
- number of species split across multiple clusters, and
- precision and recall (as described above) of the clustering results

These are the commands used to generate the evaluation results under `results/chimera_eval/`:

```bash
#constants
taxfile="data/AllSamples_2209.221216/asv_taxa.tsv"
$countsfile="data/AllSamples_2209.221216/asv_counts.tsv"
trusted_fasta="data/AllSamples_2209.221216/asv_seqs_in_finbol.fasta"

#runs
run_name="run1"
chimera_run_name="bw.strict"
clustfile="results/opticlust/AllSamples_2209.221216/$chimera_run_name/samplewise.uchime_denovo/Family/runs/$run_name/asv_reps.taxonomy.tsv"
chimera_fasta="results/chimera/AllSamples_2209.221216/filtered/$chimera_run_name/samplewise.uchime_denovo/chimeras.fasta"
nonchimera_fasta="results/chimera/AllSamples_2209.221216/filtered/$chimera_run_name/samplewise.uchime_denovo/nonchimeras.fasta"
outfile="results/chimera_eval/AllSamples_2209.221216.opticlust.$chimera_run_name.$run_name.tsv"
logfile="results/chimera_eval/AllSamples_2209.221216.opticlust.$chimera_run_name.$run_name.log"
python workflow/scripts/evaluate_chimeras.py -t $taxfile -c $clustfile --counts $countsfile --trusted_fasta $trusted_fasta \
    --chimera_fasta $chimera_fasta --nonchimera_fasta $nonchimera_fasta > $outfile 2> $logfile

run_name="run2"
chimera_run_name="bw.lenient"
clustfile="results/opticlust/AllSamples_2209.221216/$chimera_run_name/samplewise.uchime_denovo/Family/runs/$run_name/asv_reps.taxonomy.tsv"
chimera_fasta="results/chimera/AllSamples_2209.221216/filtered/$chimera_run_name/samplewise.uchime_denovo/chimeras.fasta"
nonchimera_fasta="results/chimera/AllSamples_2209.221216/filtered/$chimera_run_name/samplewise.uchime_denovo/nonchimeras.fasta"
outfile="results/chimera_eval/AllSamples_2209.221216.opticlust.$chimera_run_name.$run_name.tsv"
logfile="results/chimera_eval/AllSamples_2209.221216.opticlust.$chimera_run_name.$run_name.log"
python workflow/scripts/evaluate_chimeras.py -t $taxfile -c $clustfile --counts $countsfile --trusted_fasta $trusted_fasta \
    --chimera_fasta $chimera_fasta --nonchimera_fasta $nonchimera_fasta > $outfile 2> $logfile

run_name="run3"
chimera_run_name="sw.strict"
clustfile="results/opticlust/AllSamples_2209.221216/$chimera_run_name/samplewise.uchime_denovo/Family/runs/$run_name/asv_reps.taxonomy.tsv"
chimera_fasta="results/chimera/AllSamples_2209.221216/filtered/$chimera_run_name/samplewise.uchime_denovo/chimeras.fasta"
nonchimera_fasta="results/chimera/AllSamples_2209.221216/filtered/$chimera_run_name/samplewise.uchime_denovo/nonchimeras.fasta"
outfile="results/chimera_eval/AllSamples_2209.221216.opticlust.$chimera_run_name.$run_name.tsv"
logfile="results/chimera_eval/AllSamples_2209.221216.opticlust.$chimera_run_name.$run_name.log"
python workflow/scripts/evaluate_chimeras.py -t $taxfile -c $clustfile --counts $countsfile --trusted_fasta $trusted_fasta \
    --chimera_fasta $chimera_fasta --nonchimera_fasta $nonchimera_fasta > $outfile 2> $logfile

run_name="run4"
chimera_run_name="sw.lenient"
clustfile="results/opticlust/AllSamples_2209.221216/$chimera_run_name/samplewise.uchime_denovo/Family/runs/$run_name/asv_reps.taxonomy.tsv"
chimera_fasta="results/chimera/AllSamples_2209.221216/filtered/$chimera_run_name/samplewise.uchime_denovo/chimeras.fasta"
nonchimera_fasta="results/chimera/AllSamples_2209.221216/filtered/$chimera_run_name/samplewise.uchime_denovo/nonchimeras.fasta"
outfile="results/chimera_eval/AllSamples_2209.221216.opticlust.$chimera_run_name.$run_name.tsv"
logfile="results/chimera_eval/AllSamples_2209.221216.opticlust.$chimera_run_name.$run_name.log"
python workflow/scripts/evaluate_chimeras.py -t $taxfile -c $clustfile --counts $countsfile --trusted_fasta $trusted_fasta \
    --chimera_fasta $chimera_fasta --nonchimera_fasta $nonchimera_fasta > $outfile 2> $logfile
```

### Clustering benchmark

We benchmarked the performance of three tools for clustering ASVs into OTUs:

- SWARM (v3.1.0)
- opticlust (implemented in mothur v1.44.11)
- dbotu3 (v1.5.3)

For each tool we filtered chimeras using the strict samplewise method described above. We ran each tool with a range of parameters (see table below) and evaluated the results by calculating precision and recall values as described above.

| Run name | Swarm | Opticlust | dbOTU3 | Config file |
|----------|-------|-----------|--------|-------------|
| eval1 (default) | d=1, fastidious=True, boundary=3 | cutoff=0.03 | dist=0.1, abund=10 | `config/cluster_eval/eval1.yml` |
| eval2 | d=3 | cutoff=0.01  | dist=0.01, abund=10 | `config/cluster_eval/eval2.yml` |
| eval3 | d=5 | cutoff=0.015 | dist=0.1, abund=30 | `config/cluster_eval/eval3.yml` |
| eval4 | d=7 | cutoff=0.02 | dist=0.15, abund=10 | `config/cluster_eval/eval4.yml` |
| eval5 | d=9 | cutoff=0.025 | dist=0.01, abund=0 | `config/cluster_eval/eval5.yml` |
| eval6 | d=11 | cutoff=0.04 | dist=0.1, abund=0 | `config/cluster_eval/eval6.yml` |
| eval7 | d=13 | cutoff=0.05 | dist=0.15, abund=0 | `config/cluster_eval/eval7.yml` |
| eval8 | d=15 | cutoff=0.07 | dist=0.01, abund=20 | `config/cluster_eval/eval8.yml` |
| eval9 | d=17 | cutoff=0.1 | dist=0.1, abund=20 | `config/cluster_eval/eval9.yml` |
| eval10 | d=19 | cutoff=0.15 | dist=0.15, abund=20 | `config/cluster_eval/eval10.yml` |

These configurations were run as:

```bash
snakemake --profile slurm --configfile <config-file>
```