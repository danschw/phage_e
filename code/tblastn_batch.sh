#!/bin/bash


out_dir=/N/slate/danschw/phage_e
mkdir -p $out_dir

git_dir=~/GitHub/phage_e
bash_dir=$git_dir/code/batch_scripts
mkdir -p $bash_dir

cur_out_dir=$out_dir/$1

samples=($(find $cur_out_dir/*.fasta))




for sample in "${samples[@]}"
do

	# extract sample name
	sample_name=($(echo "${sample##*/}"))
	sample_name=($(echo "${sample_name%%.*}"))

	#output folder for blastdb
	OUT_db="$cur_out_dir/blastdb/$sample_name"
	mkdir -p $OUT_db
	
	#output folder for data
	OUT_blast="$git_dir/data/recruitment/tblastn"
	

	# write batch script
	bash_out="$bash_dir/tblastn_$1_${sample_name}.sh"

	if [ -f $bash_out ]; then
		rm -f $bash_out
	fi

    echo '#!/bin/bash' >> $bash_out
    echo '#SBATCH --mail-user=danschw@iu.edu' >> $bash_out
    echo '#SBATCH --nodes=1' >> $bash_out
    echo '#SBATCH --ntasks-per-node=1' >> $bash_out
    echo '#SBATCH --cpus-per-task=8' >> $bash_out
    echo '#SBATCH --time=1:59:00' >> $bash_out
    echo '#SBATCH --mem=20gb' >> $bash_out
    echo '#SBATCH --mail-type=FAIL,BEGIN,END' >> $bash_out
    echo "#SBATCH --job-name=${sample_name}" >> $bash_out
    echo '' >> $bash_out
    echo '##### load dependencies #####' >> $bash_out
    echo 'module load blast' >> $bash_out
    echo '' >> $bash_out
    echo '##### make blastdb #####' >> $bash_out
    echo "cd $OUT_db" >> $bash_out
    echo "makeblastdb -in $sample -dbtype 'nucl' -out $sample_name" >> $bash_out
    echo '' >> $bash_out
    echo '##### tblastn #####' >> $bash_out
    echo "mkdir -p $OUT_blast" >> $bash_out
    echo '' >> $bash_out
    echo "tblastn -query $git_dir/data/recruitment/phageE_proteins.faa -db $OUT_db/$sample_name -out $OUT_blast/${sample_name}.tsv -evalue 1e-4 -outfmt 6 -num_threads 8 " >> $bash_out
    
	sbatch $bash_out
	echo "${sample_name} submitted"
done

