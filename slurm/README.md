# Configuration profile for SLURM systems

This directory contains a configuration profile for running the workflow on HPC systems using the SLURM workload manager.

## Using the profile

To use this profile, first set `slurm_account` to your SLURM compute account id
and set `slurm_partitino` to the default SLURM partition under the
`default-resources:` section in the `config.yaml`. To see available partitions
you can run `sinfo --summarize` on the HPC system. Also check the relevant
documentation for the HPC system to see what partitions are available.

For example, if your compute account is `naiss2024-1-100` and a partition named
`core` is available to run standard jobs, the file should look like this:
    
```yaml
default-resources: 
  slurm_account: naiss2024-1-100
  slurm_partition: core
  runtime: 240
  tasks: f"{threads}"
  cpus_per_task: 1
```

Then you can run Snakemake from the repository root, adding `--profile slurm` to
the command.

Example:

```bash
snakemake --configfile <path-to-your-configfile.yml> --profile slurm <additional-arguments>
```
