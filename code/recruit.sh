#!/bin/bash



out_dir=/N/slate/danschw/phage_e
mkdir -p $out_dir

git_dir=~/GitHub/phage_e


# get reads from SRA
bash $git_dir/code/sra-batch.sh





