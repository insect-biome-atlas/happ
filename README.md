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
subdirectory under `data/` that should contain a file `asv_seqs.fasta`, a fasta
file with sequences for ASVs, and `asv_counts.tsv`, a tab-separated file with 
ASV ids as rows and sample names as columns that describes the counts of ASVs
in the different samples.

As an example, with the subdirectory `project1` under `data/` like so:

```bash
data/
├── project1
│ ├── asv_counts.tsv
│ ├── asv_seqs.fasta
```

you should set `rundir: project1`. This can be done either in a configuration
file in YAML format:

```yaml
rundir: project1
```

which you then point to with `--configfile <your-config-file>.yml` in the 
snakemake call. Or you can set it directly with `--config rundir=project1` when
you start the workflow.

## More references

- [Brandt et al 2021](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13398)