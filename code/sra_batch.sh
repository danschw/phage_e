#!/bin/bash


out_dir=/N/slate/danschw/phage_e
mkdir -p $out_dir

git_dir=~/GitHub/phage_e
bash_dir=$git_dir/code/batch_scripts
mkdir -p $bash_dir
sra_acc_dir=$git_dir/data/recruitment/sra_acc/

samples=($(find $sra_acc_dir*.txt))




for sample in "${samples[@]}"
do

	# extract sample name
	sample_name=($(echo "${sample##*/}"))
	sample_name=($(echo "${sample_name%%.*}"))

	# get acc list for sra
	mapfile -t ACC < ${sample}

	#output folder for data
	OUT_sra="$out_dir/${sample_name}"


	# write batch script
	bash_out="$bash_dir/sra_${sample_name}.sh"

	if [ -f $bash_out ]; then
		rm -f $bash_out
	fi

    echo '#!/bin/bash' >> $bash_out
    echo '#SBATCH --mail-user=danschw@iu.edu' >> $bash_out
    echo '#SBATCH --nodes=1' >> $bash_out
    echo '#SBATCH --ntasks-per-node=8' >> $bash_out
    echo '#SBATCH --time=9:59:00' >> $bash_out
    echo '#SBATCH --mem=50gb' >> $bash_out
    echo '#SBATCH --mail-type=FAIL,BEGIN,END' >> $bash_out
    echo "#SBATCH --job-name=${sample}" >> $bash_out
    echo '' >> $bash_out
    echo '##### load dependencies #####' >> $bash_out
    echo 'module load sra-toolkit' >> $bash_out
    echo '' >> $bash_out
    echo '##### run sra fastq-dump #####' >> $bash_out
    echo "mkdir -p $OUT_sra" >> $bash_out
    echo "cd ${OUT_sra}" >> $bash_out
    echo '' >> $bash_out
    echo "fastq-dump --gzip --split-files ${ACC[@]}" >> $bash_out

    sbatch $bash_out
	echo "${sample_name} submitted"
done

