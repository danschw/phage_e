#!/bin/bash
#SBATCH --mail-user=danschw@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=9:59:00
#SBATCH --mem=50gb
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=/N/u/danschw/Carbonate/GitHub/phage_e/data/recruitment/sra_acc/GOS.txt

##### load dependencies #####
module load sra-toolkit

##### run sra fastq-dump #####
mkdir -p /N/slate/danschw/phage_e/GOS
cd /N/slate/danschw/phage_e/GOS

fastq-dump --skip-technical --readids --dumpbase --fasta 60 --split-files --clip SRR9417992 SRR9417993 SRR9417994 SRR9417995 SRR9417996 SRR9417997 SRR9417998 SRR9417999 SRR9418000 SRR9418001 SRR9418002 SRR9418003 SRR9418004 SRR9418005 SRR9418006 SRR9418007 SRR9418008 SRR9418009 SRR9418010 SRR9418011 SRR9418012 SRR9418013 SRR9418014 SRR9418015 SRR9418016 SRR9418017 SRR9418018 SRR9418019 SRR9418020 SRR9418021 SRR9418022 SRR9418023 SRR9418024 SRR9418025 SRR9418026 SRR9418027 SRR9418028 SRR9418029 SRR9418030 SRR9418031 SRR9418032 SRR9418033 SRR9418034 SRR9418035 SRR9418036 SRR9418037 SRR9418038 SRR9418039 SRR9418040 SRR9418041 SRR9418042 SRR9418043 SRR9418044 SRR9418045 SRR9418046 SRR9418047 SRR9418048 SRR9418049 SRR9418050 SRR9418051 SRR9418052 SRR9418053 SRR9418054 SRR9418055 SRR9418056 SRR9418057 SRR9418058 SRR9418059 SRR9418060 SRR9418061 SRR9418062 SRR9418063 SRR9418064 SRR9418065 SRR9418066 SRR9418067 SRR9418068 SRR9418069 SRR9418070 SRR9418071 SRR9418072 SRR9418073 SRR9418074 SRR9418075 SRR9418076 SRR9418077 SRR9418078 SRR9418079 SRR9418080 SRR9418081 SRR9418082 SRR9418083 SRR9418084 SRR9418085 SRR9418086 SRR9418087 SRR9418088 SRR9418089 SRR9418090 SRR9418091 SRR066138 SRR066139 SRR192261 SRR192262 SRR192263 SRR192264 SRR192265 SRR192554 SRR192557 SRR192559 SRR192660 SRR192661 SRR192665 SRR192666 SRR192667 SRR192668 SRR192670 SRR192672 SRR192673 SRR192674 SRR192675 SRR192676 SRR341935 SRR341936 SRR341937 SRR342213 SRR9417924 SRR9417925 SRR9417926 SRR9417927 SRR9417928 SRR9417929 SRR9417930 SRR9417931 SRR9417932 SRR9417933 SRR9417934 SRR9417935 SRR9417936 SRR9417937 SRR9417938 SRR9417939 SRR9417940 SRR9417941 SRR9417942 SRR9417943 SRR9417944 SRR9417945 SRR9417946 SRR9417947 SRR9417948 SRR9417949 SRR9417950 SRR9417951 SRR9417952 SRR9417953 SRR9417954 SRR9417955 SRR9417956 SRR9417957 SRR9417958 SRR9417959 SRR9417960 SRR9417961 SRR9417962 SRR9417963 SRR9417964 SRR9417965 SRR9417966 SRR9417967 SRR9417968 SRR9417969 SRR9417970 SRR9417971 SRR9417972 SRR9417973 SRR9417974 SRR9417975 SRR9417976 SRR9417977 SRR9417978 SRR9417979 SRR9417980 SRR9417981 SRR9417982 SRR9417983 SRR9417984 SRR9417985 SRR9417986 SRR9417987 SRR9417988 SRR9417989 SRR9417990 SRR9417991 
