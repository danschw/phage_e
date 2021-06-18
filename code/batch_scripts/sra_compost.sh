#!/bin/bash
#SBATCH --mail-user=danschw@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=9:59:00
#SBATCH --mem=50gb
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=/N/u/danschw/Carbonate/GitHub/phage_e/data/recruitment/sra_acc/compost.txt

##### load dependencies #####
module load sra-toolkit

##### run sra fastq-dump #####
mkdir -p /N/slate/danschw/phage_e/compost
cd /N/slate/danschw/phage_e/compost

fastq-dump --skip-technical --readids --dumpbase --fasta 60 --split-files --clip SRR7778147 SRR7778148 SRR7778149 SRR7778150 SRR7778151 SRR7778152 SRR7778153 SRR7778154 SRR7778155 SRR7778156 SRR7778157 SRR7778158 SRR7778159 SRR7778160 SRR7778161 SRR7778162 SRR7778163 SRR7778164 
