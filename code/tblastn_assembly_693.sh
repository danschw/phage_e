#!/bin/bash
#SBATCH --mail-user=danschw@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=1:59:00
#SBATCH --mem=50gb
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=blast_693

##### load dependencies #####
module load blast

sra="SRR8088693"

git_dir="/N/u/danschw/Carbonate/GitHub/phage_e"

##### make blastdb #####
db_dir=/N/slate/danschw/phage_e/blast/$sra
mkdir -p $db_dir
cd $db_dir


makeblastdb -in "/N/slate/danschw/phage_e/lm/assembly/$sra/scaffolds.fasta" -dbtype 'nucl' -out $sra

##### tblastn #####
out_dir=$git_dir/data/assembly/tblastn
mkdir -p $out_dir

tblastn -query $git_dir/data/recruitment/phageE_proteins.faa \
-db $db_dir/$sra -out "$out_dir/$sra.tsv" -evalue 1e-4 -outfmt 6 -num_threads 8 
