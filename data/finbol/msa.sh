#!/bin/bash -l

#SBATCH -A snic2022-5-211
#SBATCH -p node -n 20
#SBATCH -t 6:00:00
#SBATCH -J MSA.finbol

module load bioinfo-tools clustalo/1.2.4

fasta="finbol.COI.fasta"
out="finbol.msa.nohmm.fasta"
hmm="PF00115.hmm"

clustalo -i $fasta -o $out --dealign -t DNA --threads 20 -l clustalo.nohmm.log -v
