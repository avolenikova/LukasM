#!/bin/bash
#usage: $ bash ThisScript.sh /PATH/to/input/dir/ /PATH/to/output/dir/
#1) PREPARE ENVIRONMENT
inputDir=$1
timeStamp=$(date +"%Y%m%d.%H%M")
outDir=${2}Interlacing-sampling_${timeStamp}
mkdir $outDir
cd $outDir
#2) CHOOSE DATA AND PRINT PBS SCRIPT
for file in "${1}"*1P.fastq.gz
do
	fileName="$(basename "$file" _1P.fastq.gz)"
	baseName="$(basename "$file")"
	scriptName=${outDir}/PBS_Interlacing-sampling_"${fileName}".sh
	(
	cat << endOfPBSscript
#!/bin/bash
#PBS -N Interl-sampl
#PBS -l select=1:ncpus=16:mem=32gb:scratch_local=250gb
#PBS -l walltime=48:00:00
#PBS -o ${outDir}/${fileName}_Interl-sampl.stdout
#PBS -e ${outDir}/${fileName}_Interl-sampl.stderr

trap 'clean_scratch' TERM EXIT

cd \${SCRATCHDIR}

#copy input files to scratchdir
scp ${file} \${SCRATCHDIR}
scp ${1}${fileName}_2P.fastq.gz \${SCRATCHDIR}

#conversion to fasta
zcat ${file} | sed -n '1~4s/^@/>/p;2~4p' > ${fileName}_1P.fasta
zcat ${fileName}_2P.fastq.gz | sed -n '1~4s/^@/>/p;2~4p' > ${fileName}_2P.fasta

#adding suffix to forward and reverse reads
for f in ${fileName}_1P.fasta;
do
awk -F " " '{if (/^>/) print \$1 "/1"; else print \$0}' $f>>${fileName}_rnmd_1P.fasta;
done

for f in ${fileName}_2P.fasta;
do
awk -F " " '{if (/^>/) print \$1 "/2"; else print \$0}' $f>>${fileName}_rnmd_2P.fasta;
done

#loading tools
module add repeatexplorerREportal
module unload python-2.7.10-gcc
module add python-2.7.5
module add debian8-compat

#interlacing
fasta_interlacer.py \
-a ${fileName}_rnmd_1P.fasta \
-b ${fileName}_rnmd_2P.fasta \
-p ${fileName}_interlaced.fasta \
-x ${fileName}_notInterlaced.fasta

#sampling (3 technical replicas)
sampleFasta.sh \
-f ${fileName}_interlaced.fasta \
-s 123654 \
-n 500000 \
-p true \
> ${fileName}_subs-500k_r1.fasta
sed 's/>/>1/g' ${fileName}_subs-500k_r1.fasta > ${fileName}_subs-500k_pref_r1.fasta

sampleFasta.sh \
-f ${fileName}_interlaced.fasta \
-s 321456 \
-n 500000 \
-p true \
> ${fileName}_subs-500k_r2.fasta
sed 's/>/>2/g' ${fileName}_subs-500k_r2.fasta > ${fileName}_subs-500k_pref_r2.fasta

sampleFasta.sh \
-f ${fileName}_interlaced.fasta \
-s 321654 \
-n 500000 \
-p true \
> ${fileName}_subs-500k_r3.fasta
sed 's/>/>3/g' ${fileName}_subs-500k_r3.fasta > ${fileName}_subs-500k_pref_r3.fasta

#concatenating
cat \
${fileName}_subs-500k_pref_r1.fasta \
${fileName}_subs-500k_pref_r2.fasta \
${fileName}_subs-500k_pref_r3.fasta \
>> ${fileName}_REx2_500k-each.fasta

#removing input files
rm \${SCRATCHDIR}/${baseName}
rm \${SCRATCHDIR}/${fileName}_2P.fastq.gz

#moving outputs to outdir
scp -r \${SCRATCHDIR}/* ${outDir} || export CLEAN_SCRATCH=true

endOfPBSscript

	) > ${scriptName}
	chmod +x ${scriptName}
#	qsub ${scriptName}
done
