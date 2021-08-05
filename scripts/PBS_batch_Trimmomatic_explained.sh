#!/bin/bash
#usage: $ bash TenhleScript.sh /absolutni/cesta/do/slozky/kde/mas/ready/ /absolutni/cesta/do/slozky/kde/chces/vysledek/ ###TOHLE MUSI BYT SPRAVNE

#1) PRIPRAVA PROSTREDI
inputDir=$1 #tohle je cesta do slozky s inputem (automaticky se doplni tou cestou co jsi zadal pri spousteni)
timeStamp=$(date +"%Y%m%d.%H%M") #datum a cas pri spusteni scriptu
outDir=${2}Trimmomatic_${timeStamp} #cesta do slozky s vysledky a nazev podslozky s touto ulohou oznacenou datem a casem
mkdir $outDir #vytvori tu slozku na vysledky ktera je specifikovala o radek vys
cd $outDir #vleze do te slozky a v ni bude tvorit nasledujici PBS scripty (pro kazdy dataset jeden)

#2) VYBER DATA A PRO KAZDY DATASET VYTVOR PBS SCRIPT
for file in "${1}"*1.fastq.gz #pro kazdy soubor v slozce kde mas ready (zadefinoval jsi nahore), ktery konci "1.fastq.gz"
do #udelej nasledujici
	fileName="$(basename "$file" _1.fastq.gz)" # jmeno datasetu je nazev souboru bez "_1.fastq.gz" (= FILENAME)
	baseName="$(basename "$file")" # nazev souboru s forward ready
	scriptName=${outDir}/PBS_Trimmomatic_"${fileName}".sh #vytvor PBS script v outputove podslozce, ktery se bude jmenovat podle nazvu datasetu
	(
	cat << endOfPBSscript #tady to zacina
#!/bin/bash
#PBS -N Trimmomatic #nazev ulohy
#PBS -l select=1:ncpus=16:mem=32gb:scratch_local=100gb #zdroje
#PBS -l walltime=48:00:00 #doba behu
#PBS -o ${outDir}/${fileName}_Trimmomatic.stdout #kam se vytvori standardni output a jak se bude jmenovat
#PBS -e ${outDir}/${fileName}_Trimmomatic.stderr #kam se vytvori standardni error a jak se bude jmenovat

trap 'clean_scratch' TERM EXIT

cd \${SCRATCHDIR} #vleze si na vzdaleny pocitac

scp ${file} \${SCRATCHDIR} #zkopiruje tam forward ready
scp ${1}${fileName}_2.fastq.gz \${SCRATCHDIR} #zkopiruje tam reverse ready
scp /storage/brno11-elixir/home/volen_a/klup/03_trimmomatic/*.fa \${SCRATCHDIR} || exit 1 #zkopiruje tam fastu s technickymi sekvencemi ###TADY MUSIS ZMENIT CESTU K TVEMU SOUBORU

module add trimmomatic-0.39 #prida modul s programem Trimmomatic

java -jar /software/trimmomatic/0.39/trimmomatic-0.39/Trimmomatic-0.39/trimmomatic-0.39.jar \ #spusti program Trimmomatic
PE -threads 15 -phred33 -trimlog ${fileName}_trimlog.txt \ # pair-end ready, paralelne 15 procesu, kodovani kvality phred33, zadefinuje jak se bude jmenovat log file
-basein ${baseName} \ #rekne programu jak se jmenuji input ready
-baseout ${fileName}_trimmed.fastq.gz \ #jak se maji jmenovat vysledky
ILLUMINACLIP:${SCRATCHDIR}Truseq_PE_frankenstein.fa:2:30:10:1:true \ #jak se maji orezat adaptory
SLIDINGWINDOW:5:20 CROP:135 HEADCROP:30 MINLEN:100 #jak se maji orezat ready aby tam nebyly zadne technicke sekvence

rm \${SCRATCHDIR}/${baseName} #smaz forwardy
rm \${SCRATCHDIR}/${fileName}_2.fastq.gz #smaz reverse

scp -r \${SCRATCHDIR}/* ${outDir} || export CLEAN_SCRATCH=true #zkopiruj vysledky a pokud neco zbylo na vzdalenem pocitaci, smaz to

endOfPBSscript #konec scriptu

	) > ${scriptName} #jmeno scriptu
	chmod +x ${scriptName} #da scriptu prava aby mohl byt spusten
	qsub ${scriptName} #posle ho do fronty uloh
done #to je vsechno :)
