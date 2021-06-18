#!/bin/bash
#SBATCH --mail-user=danschw@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=1:59:00
#SBATCH --mem=20gb
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=SRR9417935_1

##### load dependencies #####
module load blast

##### make blastdb #####
cd /N/slate/danschw/phage_e/GOS/blastdb/SRR9417935_1
makeblastdb -in /N/slate/danschw/phage_e/GOS/SRR9417935_1.fasta -dbtype 'nucl' -out SRR9417935_1

##### tblastn #####
mkdir -p /N/u/danschw/Carbonate/GitHub/phage_e/data/recruitment/tblastn

tblastn -query /N/u/danschw/Carbonate/GitHub/phage_e/data/recruitment/phageE_proteins.faa -db /N/slate/danschw/phage_e/GOS/blastdb/SRR9417935_1/SRR9417935_1 -out /N/u/danschw/Carbonate/GitHub/phage_e/data/recruitment/tblastn/SRR9417935_1.tsv -evalue 1e-4 -outfmt 6 -num_threads 8 
