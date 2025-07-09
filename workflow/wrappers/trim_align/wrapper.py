from snakemake.shell import shell

shell(
    "seqkit subseq --region {snakemake.params.codon_start}:-1 {snakemake.input.nuc} > {snakemake.output.nuc} 2>{snakemake.log}"
)
