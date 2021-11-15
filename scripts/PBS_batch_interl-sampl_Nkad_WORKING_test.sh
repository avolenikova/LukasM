#!/bin/bash
#usage: $ bash ThisScript.sh /PATH/to/input/dir/ /PATH/to/output/dir/
#1) PREPARE ENVIRONMENT
inputDir=$1 #tohle je cesta do slozky s inputy (informace z prikazove radky pri spousteni scriptu)
timeStamp=$(date +"%Y%m%d.%H%M") #dnesni datum a cas
outDir=${2}Interlacing-sampling_${timeStamp} #cesta do outputove slozky (informace z prikazove radky pri spousteni scriptu) + jmeno nove slozky s vysledky (obsahuje to datum
mkdir $outDir #vytvor slozku zadefinovanou na predchozim radku
cd $outDir #bez do teto slozky

#2) CHOOSE DATA AND PRINT PBS SCRIPT
for file in "${1}"*1P.fastq.gz #pro kazdy soubor ve slozce s inputy, ktery konci *1P.fastq.gz
do #udelej nasledujici
	fileName="$(basename "$file" _1P.fastq.gz)" #filename je nazev souboru bez koncovky _1P.fastq.gz (tj. jak se jmenuje ten dataset cely, napr. Trichonephila_10x)
	baseName="$(basename "$file")" #basename je nazev souboru (napr. Trichonephila_10x_1P.fastg.gz)
	scriptName=${outDir}/PBS_Interlacing-sampling_"${fileName}".sh #vytvor script, ktery bude v nove slozce na vysledky a bude obsahovat filename
	(
	cat << endOfPBSscript #zacni psat script, tohle je zacatek
#!/bin/bash #script je v bashi
#PBS -N Interl-sampl #takhle se jmenuje ta uloha
#PBS -l select=1:ncpus=4:mem=8gb:scratch_local=50gb #tolik chci zdroju
#PBS -l walltime=1:00:00 #takhle dlouho to pobezi
#PBS -o ${outDir}/${fileName}_Interl-sampl.stdout #sem zapis standardni output
#PBS -e ${outDir}/${fileName}_Interl-sampl.stderr #sem zapis standardni error

trap 'clean_scratch' TERM EXIT #pokud se neco stane, vycisti scratch a ukonci ulohu

cd \${SCRATCHDIR} #jdi do scratchdir

#copy input files to scratchdir 
scp ${file} \${SCRATCHDIR} #zkopiruj forward na scratchdir
scp ${1}${fileName}_2P.fastq.gz \${SCRATCHDIR} #zkopiruj reverse na scratchdir

#conversion to fasta
zcat ${file} | sed -n '1~4s/^@/>/p;2~4p' > ${fileName}_1P.fasta #cti komprimovany soubor | z kazdych ctyr radku vezmi jen prvni a druhy | zapis vysledek do filename_1P.fasta
zcat ${fileName}_2P.fastq.gz | sed -n '1~4s/^@/>/p;2~4p' > ${fileName}_2P.fasta # to same pro reverse

#adding suffix to forward and reverse reads
for f in ${fileName}_1P.fasta; #pro radek v souboru filename_1P.fasta
do #udelej nasledujici
awk -F " " '{if (/^>/) print \$1 "/1"; else print \$0}' \${f} >>${fileName}_rnmd_1P.fasta; #pokud radek zacina > (tj. je to header), pridej za nej /1 , jinak ho preskoc a vysledek zapis do filename_rnmd_1P.fasta
done #to je vsechno

for g in ${fileName}_2P.fasta; # to same pro reverse, akorat se prida k headerum /2 misto /1
do
awk -F " " '{if (/^>/) print \$1 "/2"; else print \$0}' \${g} >>${fileName}_rnmd_2P.fasta;
done

#loading tools = nacti si nasledujici moduly
module add repeatexplorerREportal
module unload python-2.7.10-gcc
module add python-2.7.5
module add debian8-compat

#interlacing = z forwardu a reversu vytvor jeden soubor, kde se budou stridat (forward a reverse z jednoho fragmentu bude vzdy za sebou)
fasta_interlacer.py \ #tohle je ten nastroj z modulu repeatexplorerREportal, ktery na to pouzijeme
-a ${fileName}_rnmd_1P.fasta \ #tohle je forward
-b ${fileName}_rnmd_2P.fasta \ #tohle je reverse
-p ${fileName}_interlaced.fasta \ #takhle se bude jmenovat vysledek: filename_interlaced.fasta
-x ${fileName}_notInterlaced.fasta #sem dej pripadne ready, ktere se nepovedlo sparovat (nemaji ten druhy z paru) - tento soubor by mel zustat prazdny

#sampling (3 technical replicas) = vytvor tri technicke pseudrepliky
sampleFasta.sh \ #tohle je ten nastroj z modulu repeatexplorerREportal, ktery na to pouzijeme
-f ${fileName}_interlaced.fasta \ #tohle je vstupni soubor, ten se stridajicimi se forwardy a reversy
-s 123 \ #seed - to je nahodne vybrane cislo, od ktereho se to odrazi pri tom pseudonahodnem vyberu 
-n 100 \ #chci 100 pseudonahodne vybranych readu => TADY UPRAVUJES POCET READU V PSEUDOREPLICE 1 :)
-p true \ #ready jsou v paru, chci oba ready z paru
> ${fileName}_subs-100_r1.fasta #zapis to do nasledujiciho souboru: filename_subs_100_r1.fasta - to cislo 100 oznacuje pocet readu, tj. taky to zmen
sed 's/>/>1/g' ${fileName}_subs-100_r1.fasta > ${fileName}_subs-100_pref_r1.fasta # na zacatek kazdeho headeru dej 1 => TOHLE TAKY ZMEN a dej tam tri nebo ctyrpismennou zkratku toho druhu a 1 (napr Trichonephila clavipes => >Tcla1 misto >1)

sampleFasta.sh \ # totez pro druhou repliku, opet potreba zmenit - stejny pocet readu jako u repliky 1, pridat jmeno druhu
-f ${fileName}_interlaced.fasta \
-s 321 \
-n 100 \ #TADY ZMENIT
-p true \
> ${fileName}_subs-100_r2.fasta
sed 's/>/>2/g' ${fileName}_subs-100_r2.fasta > ${fileName}_subs-100_pref_r2.fasta #TADY ZMENIT ('s/>/>2/g' napr. na 's/>/>Tcla2/g') - musi tam byt jmeno druhu a oznaceni pseudorepliky (2)

sampleFasta.sh \ #totez psudoreplika 3, opet nutno upravit stejne veci
-f ${fileName}_interlaced.fasta \
-s 654 \
-n 100 \ #TADY
-p true \
> ${fileName}_subs-100_r3.fasta
sed 's/>/>3/g' ${fileName}_subs-100_r3.fasta > ${fileName}_subs-100_pref_r3.fasta #TADY

#concatenating => spoj vsechny tri pseudorepliky do jednoho souboru
cat \
${fileName}_subs-100_pref_r1.fasta \
${fileName}_subs-100_pref_r2.fasta \
${fileName}_subs-100_pref_r3.fasta \
>> ${fileName}_REx2_500k-each.fasta

#removing input files => smaz puvodni soubory ze scratche
rm \${SCRATCHDIR}/${baseName}
rm \${SCRATCHDIR}/${fileName}_2P.fastq.gz

#moving outputs to outdir => zkopiruj vsechno do nove slozky na vysledky 
scp -r \${SCRATCHDIR}/* ${outDir} || export CLEAN_SCRATCH=true

endOfPBSscript #konec scriptu

	) > ${scriptName} #zapis to do souboru ktery ma zadefinovany nazev nahore (scriptName)
	chmod +x ${scriptName} #dej souboru prava aby se dal spustit
	qsub ${scriptName} #posli ulohu na vypocetni cluster
done #hotovo :)
