# ASV-clustering

## Overview

Sub-project aimed at benchmarking clustering software for ASV sequences

| Software  | Reference                                                                                      | Code                                                                    |
|-----------|------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------|
| SWARM     | [Mahé et al 2014](https://peerj.com/articles/593/)                                             | [GitHub](https://github.com/torognes/swarm)                             |
| OptiClust | [Westcott & Schloss 2017](https://journals.asm.org/doi/10.1128/mSphereDirect.00073-17)         | [GitHub](https://github.com/SchlossLab/Westcott_OptiClust_mSphere_2017) |
| dbOTU3    | [Olesen et al 2017](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0176335) | [GitHub](https://github.com/swo/dbotu3)                                 |
| LuLu      | [Guldberg et al 2017](https://www.nature.com/articles/s41467-017-01312-x)                      | [GitHub](https://github.com/tobiasgf/lulu)                              |

## Installing

Use conda to install the base environment needed to run the workflow:

```bash
conda env create -f environment.yml
```

Activate the environment:

```bash
conda activate ASV-clustering
```

## Running the workflow

There are several parameters that can be set for the different tools, but the
minimum required information needed is `rundir:` which should be the name of a 
subdirectory under `data/` that should contain:

1. a file `asv_seqs.fasta` which should be a fasta file with sequences for ASVs, 
2. a file `asv_counts.tsv` which should be a tab-separated file with ASV ids as rows and sample names as columns giving the counts of ASVs  in the different samples.
3. a file `asv_taxa.tsv` which should be a tab-separated file with taxonomic assignments for ASVs

Note that the `asv_taxa.tsv` file should have ranks in the header that match with
the configuration settings for `split_rank` and `evaluation_rank` (`Family` and `Species`
by default, respectively).

As an example, with the subdirectory `project1` under `data/` like so:

```bash
data/
├── project1
│ ├── asv_counts.tsv
│ ├── asv_seqs.fasta
│ ├── asv_taxa.tsv
```

you should set `rundir: project1`. This can be done either in a configuration
file in YAML format:

```yaml
rundir: project1
```

which you then point to with `--configfile <your-config-file>.yml` in the 
snakemake call. Or you can set it directly with `--config rundir=project1` when
you start the workflow.

### Test run
You can try the workflow out by running it on a small test dataset under `data/test/`.
To do so, simply run:

```bash
snakemake --profile test
```

## Workflow overview
The idea with this workflow is to make it easy to run OTU clustering with many 
different parameter settings then evaluate which settings you think works best
for your data. `vsearch` makes up the basis of the workflow by creating pairwise
alignments of the sequences (often the most time-consuming step) and several 
clustering runs can then be executed without having to re-create the alignments
for each run.

### 1. Splitting the input
The workflow is set up to split the input data by taxonomic rank family and 
perform clustering of sequences in each family in parallell. However, the data
can be arbitrarily split on any taxonomic rank, configurable via the `split_rank`
config parameter. The workflow will by default run on all unique taxa names at 
rank `split_rank` found in the `data/{rundir}/asv_taxa.tsv` file, and will list
the taxa in the file `data/{rundir}/{split_rank}.txt`. You can create that latter
file yourself if you want, which will override the default behaviour and cause 
the workflow to only run on taxa listed. For example, if you list (one per row)
only your families of interest in `data/{rundir}/Family.txt` then the workflow
will only run with sequences belonging to those families.
Note that if you do not specify a file with ranks, _e.g._ `Family.txt` then the
workflow will read the `asv_taxa.tsv` file, identify taxa at `{split_rank}` but 
only output those that have at least one ASV with a total count > 0 in the counts
file **and** that have a sequence in the `asv_seqs.fasta` file.

### 2. Filtering and formatting
Initially, the counts for each sequence is summed across samples and any potential 
sequence with a sum of 0 is removed. The filtering step also ensures that the 
taxonomic and sequence data lines up with no missing sequences in either.

### 3. Aligning
Filtered sequences are then aligned per `split_rank` (Family by default) using
`vsearch` like so:

```bash
vsearch --usearch_global <input> --db <input> --self --userout <output-distance> \
  -userfields query+target+id --maxaccepts 0 --maxrejects 0 --id {params.id} \
  --iddef {params.iddef}  --query_cov {params.query_cov} --threads {threads}
```

where `params.id=0.84`, `params.iddef=1`, `params.query_cov=0.9` are the default
settings (configurable via the config file).

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

```bash
results/
├── opticlust
│   └── run1
│       ├── Taxa1
│       │   └── default
│       │       └── asv_clusters.tsv
│       ├── Taxa2
│       │   └── default
│       │       └── asv_clusters.tsv
│       └── Taxa3
│           └── default
│               └── asv_clusters.tsv
├── swarm
│   └── run1
│       ├── Taxa1
│       │   └── default
│       │       └── asv_clusters.tsv
│       ├── Taxa2
│       │   └── default
│       │       └── asv_clusters.tsv
│       └── Taxa3
│           └── default
│               └── asv_clusters.tsv
└── vsearch
    └── run1
        ├── Taxa1
        │   ├── asv_seqs.dist.gz
        │   └── asv_seqs.dist.reformat.gz
        ├── Taxa2
        │   ├── asv_seqs.dist.gz
        │   └── asv_seqs.dist.reformat.gz
        └── Taxa3
            ├── asv_seqs.dist.gz
            └── asv_seqs.dist.reformat.gz
```

In the example above, the configfile used contained:
```yaml
rundir="run1"
software=["swarm", "opticlust"]
run_name="default"
```

You can see that each tool gets its own subdirectory under `results/`,
and that `rundir` (in this case "run1") in turn is created under each tool.
In addition, `run_name` is created as a subdirectory under each taxa directory. 
This way the underlying data (ASV sequences, counts and taxonomic info) found in 
the `data/<rundir>` directory can be used for any number of `run_name` settings.

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

### FINBOL

We benchmarked the workflow on the FINBOL database (see the 
[README](data/finbol/README.md)) using the following parameters:

| run name  | swarm                  | opticlust    |
|-----------|------------------------|--------------|
| default   | -d 1 --fastidious -b 3 | cutoff=0.03  |
| params1   | -d 2 -b 0              | cutoff=0.01  |
| params2   | -d 3 -b 0              | cutoff=0.015 |
| params3   | -d 4 -b 0              | cutoff=0.02  |
| params4   | -d 5 -b 0              | cutoff=0.025 |
| params5   | -d 4 -m 6 -p 3 -b 0    | cutoff=0.035 |
| params6   | -d 13 -b 0             | cutoff=0.04  |
| params7   | -d 15 -b 0             | cutoff=0.05  |
| params8   | -d 17 -b 0             | cutoff=0.06  |
| params9   | -d 20 -b 0             | cutoff=0.07  |
| params10  | -d 22 -b 0             | cutoff=0.08  |
| params11  | -d 23 -b 0             | cutoff=0.09  |
| params12  | -d 25 -b 0             | cutoff=0.1   |

Our analysis showed that `params4` with opticlust (`cutoff=0.025`) gave the 
highest precision/recall. However, when we used this setting on real-world data
for >630,000 ASVs in >5000 samples the results showed a very low recall and 
a much greater number of clusters than anticipated. Looking more closely at the
Lepidoptera order there were ~4,300 clusters but only ~1,400 unique species
in the data. We hypothesized that chimeric ASV sequences were causing problems
for the clustering algorithms, something we did not observe in the benchmarking 
because it was done on presumably high-quality non-chimeric sequences.


## More references

- [Brandt et al 2021](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13398)