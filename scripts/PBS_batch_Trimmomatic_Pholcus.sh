#!/bin/bash
#usage: $ bash ThisScript.sh /PATH/to/input/dir/ /PATH/to/output/dir/ ###TOHLE MUSI BYT SPRAVNE
#1) PREPARE ENVIRONMENT
inputDir=$1
timeStamp=$(date +"%Y%m%d.%H%M")
outDir=${2}Trimmomatic_${timeStamp}
mkdir $outDir
cd $outDir
#2) CHOOSE DATA AND PRINT PBS SCRIPT
for file in "${1}"*1.fastq.gz
do
	fileName="$(basename "$file" _1.fastq.gz)"
	baseName="$(basename "$file")"
	scriptName=${outDir}/PBS_Trimmomatic_"${fileName}".sh
	(
	cat << endOfPBSscript
#!/bin/bash
#PBS -N Trimmomatic
#PBS -l select=1:ncpus=16:mem=32gb:scratch_local=100gb
#PBS -l walltime=48:00:00
#PBS -o ${outDir}/${fileName}_Trimmomatic.stdout
#PBS -e ${outDir}/${fileName}_Trimmomatic.stderr

trap 'clean_scratch' TERM EXIT

cd \${SCRATCHDIR}

scp ${file} \${SCRATCHDIR}
scp ${1}${fileName}_2.fastq.gz \${SCRATCHDIR}
scp /storage/brno11-elixir/home/volen_a/klup/03_trimmomatic/*.fa \${SCRATCHDIR} || exit 1 ###TADY UPRAV CESTU

module add trimmomatic-0.39

java -jar /software/trimmomatic/0.39/trimmomatic-0.39/Trimmomatic-0.39/trimmomatic-0.39.jar \
PE -threads 15 -phred33 -trimlog ${fileName}_trimlog.txt \
-basein ${baseName} \
-baseout ${fileName}_trimmed.fastq.gz \
ILLUMINACLIP:${SCRATCHDIR}Truseq_PE_frankenstein.fa:2:30:10:1:true \
SLIDINGWINDOW:5:20 CROP:135 HEADCROP:30 MINLEN:100

rm \${SCRATCHDIR}/${baseName}
rm \${SCRATCHDIR}/${fileName}_2.fastq.gz

scp -r \${SCRATCHDIR}/* ${outDir} || export CLEAN_SCRATCH=true

endOfPBSscript

	) > ${scriptName}
	chmod +x ${scriptName}
	qsub ${scriptName}
done
