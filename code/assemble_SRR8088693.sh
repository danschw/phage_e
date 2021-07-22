#!/bin/bash
#SBATCH --mail-user=danschw@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --time=11:59:00
#SBATCH --mem=250gb
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=SRR8088693

##### load dependencies #####
module load sra-toolkit
module load spades # v3.14.1

sra="SRR8088693"

##### run sra fastq-dump #####
mkdir -p /N/slate/danschw/phage_e/$sra
cd /N/slate/danschw/phage_e/$sra

fastq-dump --skip-technical --readids --dumpbase --split-files --clip $sra

reads=($(find *.fastq))

out_dir="/N/slate/danschw/phage_e/assembly/$sra"
mkdir -p $out_dir

spades.py --meta -t 24 -m 250 --checkpoints 'last' \
-1 ${reads[0]} \
-2 ${reads[1]} \
--only-assembler \
--phred-offset 33 \
-o $out_dir 

#cleanup
rm -r /N/slate/danschw/phage_e/$sra
       



