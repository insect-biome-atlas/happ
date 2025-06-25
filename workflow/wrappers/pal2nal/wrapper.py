from snakemake.shell import shell

shell(
    "pal2nal.pl {snakemake.input.pep} {snakemake.input.nuc} -output fasta -codontable {snakemake.params.codon_table} > {snakemake.output} 2>{snakemake.log}"
)
