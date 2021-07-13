#!/bin/bash



out_dir=/N/slate/danschw/phage_e
mkdir -p $out_dir

git_dir=~/GitHub/phage_e

bash_dir=$git_dir/code/batch_scripts
mkdir -p $bash_dir
sra_acc_dir=$git_dir/data/recruitment/sra_acc/

sample=($(find $sra_acc_dir*Kinneret.txt))


# extract sample name
sample_name=($(echo "${sample##*/}"))
sample_name=($(echo "${sample_name%%.*}"))

# get acc list for sra
mapfile -t ACC < ${sample}



for sample in "${ACC[@]}"
do

	cur_out_dir=$out_dir/$sample	

	#output folder for sra data
	OUT_sra="$cur_out_dir"
	mkdir -p $OUT_sra
	
	#output folder for blastdb
	OUT_db="$cur_out_dir/blastdb"
	
	#output folder for blast results data
	OUT_blast="$git_dir/data/recruitment/tblastn-$sample_name"

	
	# write batch script
	bash_out="$bash_dir/sra_${sample_name}_$sample.sh"

	if [ -f $bash_out ]; then
		rm -f $bash_out
	fi

    echo '#!/bin/bash' >> $bash_out
    echo '#SBATCH --mail-user=danschw@iu.edu' >> $bash_out
    echo '#SBATCH --nodes=1' >> $bash_out
    echo '#SBATCH --ntasks-per-node=1' >> $bash_out
    echo '#SBATCH --cpus-per-task=8' >> $bash_out
    echo '#SBATCH --time=4:59:00' >> $bash_out
    echo '#SBATCH --mem=20gb' >> $bash_out
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
    echo "fastq-dump --skip-technical --readids --dumpbase --fasta 60 --split-files --clip ${sample}" >> $bash_out
       # https://edwards.sdsu.edu/research/fastq-dump/
    echo 'fastas=($(find ~+ ./*.fasta))'>> $bash_out #full path with ~+

	# blast
    echo '' >> $bash_out
    echo '' >> $bash_out
    echo '# Recruitment with tblastn' >> $bash_out
    echo '##### load dependencies #####' >> $bash_out
    echo 'module load blast' >> $bash_out
    echo '' >> $bash_out
    echo '##### make blastdb #####' >> $bash_out
    echo "mkdir -p $OUT_db" >> $bash_out
    echo "cd $OUT_db" >> $bash_out
    echo 'cat ${fastas[@]} | makeblastdb -dbtype ''nucl'' -out' "$sample" ' -title' "$sample">> $bash_out
    echo '' >> $bash_out
    echo '##### tblastn #####' >> $bash_out
    echo "mkdir -p $OUT_blast" >> $bash_out
    echo '' >> $bash_out
    echo "tblastn -query $git_dir/data/recruitment/phageE_proteins.faa -db $OUT_db/$sample -out $OUT_blast/${sample_name}_$sample.tsv -evalue 1e-4 -outfmt 6 -num_threads 8 " >> $bash_out
    echo '' >> $bash_out
    echo '#cleanup' >> $bash_out
    echo "rm -r $OUT_sra" >> $bash_out



	sbatch $bash_out
done
    
########################



