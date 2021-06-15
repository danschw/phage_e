#!/bin/bash
#SBATCH --mail-user=danschw@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=9:59:00
#SBATCH --mem=50gb
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=/N/u/danschw/Carbonate/GitHub/phage_e/data/recruitment/sra_acc/subsurface.txt

##### load dependencies #####
module load sra-toolkit

##### run sra fastq-dump #####
mkdir -p /N/slate/danschw/phage_e/subsurface
cd /N/slate/danschw/phage_e/subsurface

fastq-dump --gzip --split-files SRR5716301 SRR5716302 SRR5716303 SRR5716304 SRR5716305 SRR5716306 SRR5716307 SRR5716308 SRR5716309 SRR5716310 SRR5716311 SRR5716312 SRR5716313 SRR5716314 
