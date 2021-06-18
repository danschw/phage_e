#!/bin/bash
#SBATCH --mail-user=danschw@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=9:59:00
#SBATCH --mem=50gb
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=/N/u/danschw/Carbonate/GitHub/phage_e/data/recruitment/sra_acc/groundwaters.txt

##### load dependencies #####
module load sra-toolkit

##### run sra fastq-dump #####
mkdir -p /N/slate/danschw/phage_e/groundwaters
cd /N/slate/danschw/phage_e/groundwaters

fastq-dump --skip-technical --readids --dumpbase --fasta 60 --split-files --clip SRR10598216 SRR10598217 SRR10598218 SRR10598219 SRR10598220 SRR10598221 SRR10598222 SRR10598223 SRR10598224 
