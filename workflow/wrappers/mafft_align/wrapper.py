from snakemake.shell import shell

shell(
    "mafft --auto --thread {snakemake.threads} {snakemake.input} > {snakemake.output} 2>{snakemake.log}"
)
