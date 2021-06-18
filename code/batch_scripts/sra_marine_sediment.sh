#!/bin/bash
#SBATCH --mail-user=danschw@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=9:59:00
#SBATCH --mem=50gb
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=/N/u/danschw/Carbonate/GitHub/phage_e/data/recruitment/sra_acc/marine_sediment.txt

##### load dependencies #####
module load sra-toolkit

##### run sra fastq-dump #####
mkdir -p /N/slate/danschw/phage_e/marine_sediment
cd /N/slate/danschw/phage_e/marine_sediment

fastq-dump --skip-technical --readids --dumpbase --fasta 60 --split-files --clip SRR1555743 SRR1555744 SRR1555748 SRR1555750 
