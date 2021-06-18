#!/bin/bash



out_dir=/N/slate/danschw/phage_e
mkdir -p $out_dir

git_dir=~/GitHub/phage_e


# get reads from SRA
bash $git_dir/code/sra-batch.sh

#--------------

# make blastdb and run tblastn
	# code accepts argument of directory name containing fasta files obtained in previous step
bash $git_dir/code/tblastn_batch.sh marine_sediment
bash $git_dir/code/tblastn_batch.sh subsurface
bash $git_dir/code/tblastn_batch.sh GOS
bash $git_dir/code/tblastn_batch.sh compost
bash $git_dir/code/tblastn_batch.sh groundwaters
#--------------

tblastn -query $git_dir/data/recruitment/phageE_proteins.faa -db $out_dir/test/SRR7778149 -out $out_dir/test/outputfile.tsv -evalue 1e-4 -outfmt 6 -num_threads 8 
## Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score

cut -f 2 $out_dir/test/outputfile.tsv > $out_dir/test/hit.list

#Extract sequences with names in file name.lst, one sequence name per line:
gunzip -c $out_dir/test/SRR7778149*.gz | \
grep -w -A 1 -f  $out_dir/test/hit.list --no-group-separator \
> $out_dir/test/hits.fna
