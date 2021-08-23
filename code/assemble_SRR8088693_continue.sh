#!/bin/bash
#SBATCH --mail-user=danschw@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=11:59:00
#SBATCH --mem=500gb
#SBATCH --partition=largememory
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=lm693

##### load dependencies #####
#module load sra-toolkit
module load spades # v3.14.1

sra="SRR8088693"

##### run sra fastq-dump #####
mkdir -p /N/slate/danschw/phage_e/lm/$sra
cd /N/slate/danschw/phage_e/lm/$sra

#fastq-dump --skip-technical --readids --dumpbase --split-files --clip $sra

reads=($(find *.fastq))

out_dir="/N/slate/danschw/phage_e/lm/assembly/$sra"
mkdir -p $out_dir

spades.py -t 16 -m 500 --checkpoints 'last' --restart-from 'last' \
-o $out_dir 

#cleanup
#rm -r /N/slate/danschw/phage_e/lm/$sra
       



