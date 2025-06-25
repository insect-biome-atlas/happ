from snakemake.shell import shell
import os
import shutil

tmpdir = os.path.expandvars(snakemake.params.tmpdir)

os.makedirs(tmpdir, exist_ok=True)
with open(f"{tmpdir}/cluster_reps.fasta", "w") as fhout, open(
    snakemake.input[0], "r"
) as fhin:
    for line in fhin:
        line = line.rstrip()
        if line.startswith(">"):
            line = ">" + line.split(" ")[1]
        fhout.write(line + "\n")
shell(
    "vsearch --usearch_global {tmpdir}/cluster_reps.fasta --db {tmpdir}/cluster_reps.fasta --self --id .84 --iddef 1 "
    "--userout {snakemake.output} -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits {snakemake.params.maxhits} --threads {snakemake.threads} > {snakemake.log} 2>&1"
)
shutil.rmtree(tmpdir)
