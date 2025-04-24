#!/bin/bash

#export PATH="~/NANOFMDV-PIPELINE/:${PATH}"
#export NANOFMDV_ENV=~/NANOFMDV-PIPELINE

echo ""
echo "nanofmdv_single_barcode alignment script started"
echo ""

# n = nanofmdv path (or set as environment variable)
while getopts "n:r:t:q:p:m:d:x" TEST; do
	case $TEST in

	n) nanoPath=$OPTARG
	;;
	r) refPath=$OPTARG
	;;
	t) threads=$OPTARG
	;;
	q) fqName=$OPTARG
	;;
	p) primersFile=$OPTARG
	;;
	m) medModel=$OPTARG
	;;
	d) depArg=$OPTARG
	;;
	x) OPT_X=1
	;;
	esac
done


echo "nanofdmv - checking input options and arguments"

if [[ -n "${nanoPath}" ]]; then
	# remove trailing / if exists
	if [[ "$nanoPath" =~ '/'$ ]]; then 
		nanoPath=$(echo ${nanoPath} | sed 's/\/$//')	
	fi
	
	echo "option - nanofdmv path set via -n option to = ${nanoPath}"

	if [ ! -d "$nanoPath" ]; then
    	echo "ERROR - the nanofmdv directory supplied via the -n option does not exist: ${nanoPath}"
        echo "EXITING"
        exit 1
    fi
else
    if [[ -n "${NANOFMDV_ENV}" ]]; then
		nanoPath="${NANOFMDV_ENV}"
		nanoPath=$(echo ${nanoPath} | sed 's/\/$//')

		echo "option - NANOFMDV_ENV path is set - nanofmdv path = ${nanoPath}"

		if [ ! -d "$nanoPath" ]; then
	    	echo "ERROR - the nanofmdv directory supplied via the NANOFMDV_ENV does not exist: ${nanoPath}"
            echo "EXITING"
            exit 1
        fi
	else
		echo "ERROR - there is no NANOFMDV_ENV environment variable set or path supplied"
		echo "set a NANOFMDV_ENV variable or use the -n option to specify the location of the nanofmdv scripts folder"
		echo "EXITING"
		exit 1
	fi
fi



if [[ -n "${refPath}" ]]; then
	refPath=$(echo ${refPath} | sed 's/\/$//')                
        
    echo "option - nanofmdv-ref path set via -r option to = ${refPath}"

	if [ ! -d "$refPath" ]; then
        echo "ERROR - the nanofmdv-ref directory supplied via the -r option does not exist: ${refPath}"
		echo "EXITING"
		exit 1
	fi
else	
	refPath="${nanoPath}/Refs"
	echo "option - nanofmdv-ref path set to default nanofmdv-path/Refs = ${refPath}"

	if [ ! -d "$refPath" ]; then
        echo "ERROR - the nanofmdv-ref directory (derived from the nanofmdv path/env) does not exist: ${refPath}"
		echo "EXITING"
		exit 1
    fi
fi



if [ -n "$OPT_X" ]; then
	echo "option - not running medaka to polish consensus as -x option is set"
else
if [ -n "$medModel" ]; then
		echo "option - medaka variant model set via -m option to = ${medModel}"
		# unable to check currently if correct - can we check with medaka?
	else
		echo "WARNING - the medaka variant model was not specified - defaulting to r1041_e82_400bps_sup_variant_v4.3.0"
		echo "WARNING - you really should specify the correct medaka variant model - polished consensus may be incorrect"
		medModel="r1041_e82_400bps_sup_variant_v4.3.0"
	fi
fi



if [ -n "$primersFile" ]; then
	if [ ! -f "$primersFile" ]; then
		echo "ERROR - the primers file supplied with the -p option does not exist: ${primersFile}"
		echo "EXITING"
		exit 1
	else
		pNum=$(grep -c "^>" ${primersFile})
		echo "option - primers file set via -p option: ${primersFile} - number of sequences: ${pNum}"
	fi
else
	echo "WARNING - no primers file was supplied with the -p option"
	primersFile=${nanoPath}/primers.fasta
	echo "WARNING - defaulting to standard L primers: ${nanoPath}/primers.fasta"
		
	#echo "ERROR - a primers file must be supplied with the -p option"
	#echo "EXITING"
	#exit 1
fi



if [ -n "$fqName" ]; then
	if [ ! -f "$fqName" ]; then
		echo "ERROR - the fastq file supplied with the -q option does not exist: ${fqName}"
		echo "EXITING"
		exit 1
	else
		echo "option - input fastq file set via -q option: ${fqName}"
		fqExt="${fqName##*.}"
		
		if [ "$fqExt" != fastq ] && [ "$fqExt" != fq ]; then
			echo "ERROR - the supplied fastq file does not end with .fastq or .fq"
			echo "EXITING"
			exit 1
		fi
		
		sampleFName=$(basename "${fqName}")
		sPath=$(dirname "${fqName}")
		#sName="${fqName%.f*q}"
		sName="${sampleFName%.f*q}"

		#echo "input fastq = $fqName"
		#echo "sampleFileName (no path) = ${sampleFName}"
		#echo "sampleExtension = ${fqExt}"
		#echo "sampleDir = ${sPath}"
		#echo "sample name = $sName"
	fi
else
		echo "ERROR - a fastq file must be supplied with the -q option"
		echo "EXITING"
		exit 1
fi


if [ -n "$threads" ]; then
	echo "option - threads set via the -t option to: ${threads}"
	#TODO - SANITY CHECK THREADS VALUE
else
	threads=8
	echo "option - threads -t option not supplied so defaulting to threads = ${threads}"
fi



if [ -n "$depArg" ]; then
	echo "option - minimum depth for consensus set via tha (-d) option to: ${depArg}"
	#TODO - SANITY CHECK VALUE
else
	depArg=20
	echo "option - minimum depth for cosensus (-d) option not supplied so defaulting to = ${depArg}"
fi


echo "nanofdmv - finished checking input options and arguments"

echo ""
echo "nanofdmv - pipeline started"
echo "nanofdmv - input FASTQ file = ${fqName}"
echo "nanofdmv - output filenames will all begin with = ${sName}"
echo "nanofdmv - output folder location = ${sPath}/"


rm -f ${sPath}/${sName}_vp1_idx.txt
rm -f ${sPath}/${sName}_complete_idx.txt
rm -f ${sPath}/${sName}_fmdv_log.txt
rm -f ${sPath}/${sName}_nanofmdv.log

touch ${sPath}/${sName}_vp1_idx.txt
touch ${sPath}/${sName}_complete_idx.txt
touch ${sPath}/${sName}_fmdv_log.txt
touch ${sPath}/${sName}_nanofmdv.log

echo "nanofmdv - input fastq = ${fqName}"  >> ${sPath}/${sName}_fmdv_log.txt

# align to vp1 only
echo "nanofdmv - aligning reads to a collection of fmdv VP1 sequences"

minimap2 -t ${threads} -x map-ont ${refPath}/fmdv_vp1.fasta -a ${fqName} > ${sPath}/${sName}_vp1.sam 2>>${sPath}/${sName}_nanofmdv.log
samtools sort -@${threads} ${sPath}/${sName}_vp1.sam -o ${sPath}/${sName}_vp1.bam 2>>${sPath}/${sName}_nanofmdv.log
rm ${sPath}/${sName}_vp1.sam
samtools index ${sPath}/${sName}_vp1.bam
samtools idxstats ${sPath}/${sName}_vp1.bam | sort -n -k3,3 > ${sPath}/${sName}_vp1_idx.txt

topvp1=$( tail -n1 ${sPath}/${sName}_vp1_idx.txt | cut -f1 )
vp1c=$( tail -n1 ${sPath}/${sName}_vp1_idx.txt | cut -f3 )
vp1Count=$(($vp1c))
cutoff=0

if [[ $vp1Count -gt $cutoff ]]; then
	echo "nanofdmv - top VP1 read count = ${vp1Count}, top VP1 ref = ${topvp1}"
	echo "nanofdmv - top VP1 read count = ${vp1Count}, top VP1 ref = ${topvp1}" >> ${sPath}/${sName}_fmdv_log.txt
else
	echo "WARNING"
	echo "nanofdmv - WARNING - no reads mapped to any VP1 seqs"
	echo "WARNING"
	echo "nanofdmv - WARNING - no reads mapped to any VP1 seqs" >> ${sPath}/${sName}_fmdv_log.txt
fi

# align to complete genome
echo "nanofdmv - aligning reads to a collection of fmdv complete genome sequences"

minimap2 -t ${threads} -x map-ont ${refPath}/fmdv.fasta -a ${fqName} > ${sPath}/${sName}_complete.sam 2>>${sPath}/${sName}_nanofmdv.log
samtools sort -@${threads} ${sPath}/${sName}_complete.sam -o ${sPath}/${sName}_complete.bam 2>>${sPath}/${sName}_nanofmdv.log
rm ${sPath}/${sName}_complete.sam
samtools index ${sPath}/${sName}_complete.bam
samtools idxstats ${sPath}/${sName}_complete.bam | sort -n -k3,3 > ${sPath}/${sName}_complete_idx.txt
topgen=$( tail -n1 ${sPath}/${sName}_complete_idx.txt | cut -f1 )
topgenc=$( tail -n1 ${sPath}/${sName}_complete_idx.txt | cut -f3 )
topgenCount=$(($topgenc))
cutoff=0

if [[ $topgenCount -gt $cutoff ]]; then
	if [[ ! $vp1Count -gt $cutoff ]]; then
		echo "nanofdmv - top complete genome read count = ${topgenCount}, top genome ref = ${topgen}"
		echo "nanofdmv - top complete genome read count = ${topgenCount}, top genome ref = ${topgen}" >> ${sPath}/${sName}_fmdv_log.txt
	fi
else
	echo "WARNING"
	echo "nanofmdv - WARNING - no reads mapped to any FMDV genomes" 
	echo "nanofmdv - WARNING - no reads mapped to any FMDV genomes" >> ${sPath}/${sName}_fmdv_log.txt
	echo "WARNING"
fi

echo "nanofdmv - selecting fmdv reference sequence to use"

topvp1=$( tail -n1 ${sPath}/${sName}_vp1_idx.txt | cut -f1 )
vp1c=$( tail -n1 ${sPath}/${sName}_vp1_idx.txt | cut -f3 )
vp1Count=$(($vp1c))
cutoff=0

if [[ $vp1Count -gt $cutoff ]]; then
	ref="$(tail -n1 ${sPath}/${sName}_vp1_idx.txt | cut -f1)" 
	echo "nanofmdv - SELECTED - top ref from VP1 stats = ${ref}"
	echo "nanofmdv - SELECTED - top ref from VP1 stats = ${ref}" >> ${sPath}/${sName}_fmdv_log.txt
	grep -A1 ${ref} ${refPath}/fmdv.fasta > ${sPath}/top_ref.fasta
else
	echo "nanofmdv - VP1 count was 0 - selecting reference via full genome mappings"

	topgen=$( tail -n1 ${sPath}/${sName}_complete_idx.txt | cut -f1 )
	topgenc=$( tail -n1 ${sPath}/${sName}_complete_idx.txt | cut -f3 )
	topgenCount=$(($topgenc))
	cutoff=0
        
	if [[ $topgenCount -gt $cutoff ]]; then
		ref="$(tail -n1 ${sPath}/${sName}_complete_idx.txt | cut -f1)"
		echo "nanofmdv - SELECTED - top ref via complete genome = ${ref}"
		echo "nanofmdv - SELECTED - top ref via complete genome = ${ref}" >> ${sPath}/${sName}_fmdv_log.txt
		grep -A1 ${ref} ${refPath}/fmdv.fasta > ${sPath}/top_ref.fasta
	else
		#echo "nanofmdv - WARNING - no reads mapped to any FMDV genome - defaulting to O_KC503937_2010_Andong as selected genome"
		#echo "nanofmdv - WARNING - no reads mapped to any FMDV genome - defaulting to O_KC503937_2010_Andong as selected genome" >> ${sPath}/${sName}_fmdv_log.txt
		#ref="O_KC503937_2010_Andong"
		#grep -A1 ${ref} ${refPath}/fmdv.fasta > ${sPath}/top_ref.fasta
		
		echo "nanofmdv - WARNING - no reads mapped to any FMDV genome - EXITING"
		echo "nanofmdv - WARNING - no reads mapped to any FMDV genome - EXITING" >> ${sPath}/${sName}_fmdv_log.txt
		exit 1
    fi
    
    if [[ $topgenCount -lt $depArg ]]; then
    	echo "WARNING"
		echo "nanofdmv - WARNING - the number of mapped reads (${topgenCount}) against the top FMDV genome is less than the consensus coverage threshold (${depArg}) itself- results (if any) will likely be INPRECISE"
  		echo "WARNING"
    fi
fi

rm -f ${sPath}/${sName}_top_idx.txt
touch ${sPath}/${sName}_top_idx.txt

# align the reads to the top ref seq
echo "nanofmdv - aligning reads to the selected fmdv reference genome"
minimap2 ${sPath}/top_ref.fasta -t ${threads} -a ${fqName} > ${sPath}/${sName}_top.sam 2>>${sPath}/${sName}_nanofmdv.log
samtools view -F2308 -h ${sPath}/${sName}_top.sam > ${sPath}/temp.sam
mv ${sPath}/temp.sam ${sPath}/${sName}_top.sam
samtools sort -@${threads} ${sPath}/${sName}_top.sam -o ${sPath}/${sName}_top.bam 2>>${sPath}/${sName}_nanofmdv.log
samtools index ${sPath}/${sName}_top.bam
rm -f ${sPath}/${sName}_top.sam
samtools idxstats ${sPath}/${sName}_top.bam > ${sPath}/${sName}_top_idx.txt

# use goprime to calculate primer binding and create a bed file
echo "nanofmdv - evaluating primer binding sites"
python ${nanoPath}/goprime_comp.py ${primersFile} ${sPath}/top_ref.fasta >> ${sPath}/${sName}_nanofmdv.log
python ${nanoPath}/goprime_bed.py ${sPath}/top_ref_bind.csv >> ${sPath}/${sName}_nanofmdv.log
sed "s!fmdvref!${ref}!g" ${sPath}/top_ref_bind.bed > ${sPath}/temp.bed
mv ${sPath}/temp.bed ${sPath}/top_ref_bind.bed
        
# use samtools ampliconclip to clip reads
echo "nanofmdv - clipping primers from BAM file"
samtools ampliconclip --tolerance 50 --hard-clip --both-ends --no-excluded --filter-len 50 -O SAM -b ${sPath}/top_ref_bind.bed ${sPath}/${sName}_top.bam -o ${sPath}/${sName}_clip.sam 2>>${sPath}/${sName}_nanofmdv.log
samtools view -h -F2308 ${sPath}/${sName}_clip.sam | samtools sort -@${threads} -o ${sPath}/${sName}_clip.bam 2>>${sPath}/${sName}_nanofmdv.log
samtools index ${sPath}/${sName}_clip.bam
rm -f ${sPath}/${sName}_clip.sam
samtools idxstats ${sPath}/${sName}_clip.bam > ${sPath}/${sName}_clip_idx.txt

# extract the clipped fastq reads from the BAM
echo "nanofmdv - creating primer clipped FASTQ read files" 
samtools fastq -F4 ${sPath}/${sName}_clip.bam > ${sPath}/${sName}_clip.fastq 2>>${sPath}/${sName}_nanofmdv.log

# draft consensus
echo "nanofmdv - creating mpileup for draft consensus sequence - this is normally the longest step"
samtools mpileup -aa -d0 -q 0 -Q 0 -C 0 -f ${sPath}/top_ref.fasta ${sPath}/${sName}_clip.bam > ${sPath}/${sName}_clip_mpile.txt 2>>${sPath}/${sName}_nanofmdv.log

echo "nanofmdv - creating draft consensus sequence"
#Sept2024 - cov 20 -> 10
java -jar ${nanoPath}/VSENSUS.jar ${sPath}/${sName}_clip_mpile.txt -c=${depArg} > ${sPath}/${sName}_clip_mpile_vsensus_log.txt
echo "" >> ${sPath}/${sName}_clip_mpile_vsensus.fasta
mv ${sPath}/${sName}_clip_mpile_vsensus.fasta ${sPath}/${sName}_fmdv_preconsensus.fasta
sed "s/^>.*/>${sName}_preconsensus/g" ${sPath}/${sName}_fmdv_preconsensus.fasta > ${sPath}/temp.fasta
mv ${sPath}/temp.fasta ${sPath}/${sName}_fmdv_preconsensus.fasta

preconNonN=$(grep -v ">" ${sPath}/${sName}_fmdv_preconsensus.fasta | sed 's/N//g' | wc -c | sed 's/ //g')
preconL=$(grep -v ">" ${sPath}/${sName}_fmdv_preconsensus.fasta | wc -c | sed 's/ //g')
echo "nanofmdv - preconsensus length = ${preconL} with non-N base length = ${preconNonN}"

if [ ! -n "$OPT_X" ]; then
	echo "nanofmdv - creating medaka polished consensus sequence"
	#medaka_consensus -i ${sPath}/${sName}_clip.fastq -d ${sPath}/${sName}"_fmdv_preconsensus.fasta" -o ${sPath}/medaka -t ${threads} -m ${medModel} >> ${sPath}/${sName}_nanofmdv.log 2>&1
	#cp ${sPath}/medaka/consensus.fasta ${sPath}/${sName}_fmdv_consensus.fasta
	
	medaka_variant -i ${sPath}/${sName}_clip.fastq -o ${sPath}/medaka-variant -m ${medModel} -r ${sPath}/top_ref.fasta -f -x >> ${sPath}/${sName}_nanofmdv.log 2>&1
	medaka sequence --min_depth ${depArg} --fill_char N ${sPath}/medaka-variant/consensus_probs.hdf ${sPath}/top_ref.fasta ${sPath}/medaka-variant/${sName}_fmdv_consensus.fasta >> ${sPath}/${sName}_nanofmdv.log 2>&1
	sed "s/^>.*/>${sName}_medaka_consensus/g" ${sPath}/medaka-variant/${sName}_fmdv_consensus.fasta > ${sPath}/${sName}_fmdv_consensus.fasta
	conNonN=$(grep -v ">" ${sPath}/${sName}_fmdv_consensus.fasta | sed 's/N//g' | wc -c | sed 's/ //g')
	conL=$(grep -v ">" ${sPath}/${sName}_fmdv_consensus.fasta | wc -c | sed 's/ //g')
fi

echo "nanofmdv - blasting consensus sequence"

if [ ! -n "$OPT_X" ]; then
	blastn -db ${refPath}/BLAST/fmdv-gbbase.fasta -query ${sPath}/${sName}_fmdv_consensus.fasta -out ${sPath}/${sName}_fmdv_consensus_blast.txt -outfmt 6
fi

blastn -db ${refPath}/BLAST/fmdv-gbbase.fasta -query ${sPath}/${sName}_fmdv_preconsensus.fasta -out ${sPath}/${sName}_fmdv_preconsensus_blast.txt -outfmt 6


# stats
echo "nanofmdv - generating coverage stats"
	
echo "Sample,TopConsensusBLAST,Reads,FMDVreads,TopFMDV,VP1reads,TopVP1,MappedReads,SelectedRef,RefLength,AverageDepth,Gen%Cov=0,Gen%Cov>=1,Gen%Cov>=3,Gen%Cov>=5,Gen%Cov>=10,Gen%Cov>=20,Gen%Cov>=40,Gen%Cov>=50,Gen%Cov>=100,Gen%Cov>=500,Gen%Cov>=1000,Gen%Cov>=5000,Gen%Cov>=10000,Gen%Cov>=50000,Gen%Cov>=100000,MappedReads,MappedAlignments,SecondaryAlignments,SupplementaryAlignments,UnmappedReads,BLAST-SubjectName,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore" > ${sPath}/${sName}_fmdv_stats.csv

reads=$(python ${nanoPath}/read_count.py ${fqName})
#reads=$(expr `(wc -l ${fqName}|cut -f1 -d " ")` / 4)

python ${nanoPath}/depth_charge.py ${sPath}/${sName}_clip.bam full > /dev/null
# below created by depth charge
mv ${sPath}/${sName}_clip_depth.txt ${sPath}/${sName}_fmdv_depth.txt
samtools depth -aa ${sPath}/${sName}_top.bam > ${sPath}/${sName}_top_depth.txt

echo "nanofmdv - raw read count = ${reads}" >> ${sPath}/${sName}_fmdv_log.txt

topvp1=$(cat ${sPath}/${sName}_top_idx.txt | head -n1 | cut -f1)
topgen=$(tail -n1 ${sPath}/${sName}_complete_idx.txt | cut -f1)  

map=$(samtools view -c -F2308 ${sPath}/${sName}_top.bam)
mapvp1=$(samtools view -c -F2308 ${sPath}/${sName}_vp1.bam)
mapall=$(samtools view -c -F2308 ${sPath}/${sName}_complete.bam)

echo "nanofmdv - total mapped FMDV reads (all genomes) = ${mapall}" >> ${sPath}/${sName}_fmdv_log.txt
echo "nanofmdv - mapped FMDV reads (top genome) = ${map}" >> ${sPath}/${sName}_fmdv_log.txt
mapPer=$(python -c "print( ${map} / float(${reads}) * 100)" )
echo "nanofmdv - mapped reads percentage = ${mapPer} %" >> ${sPath}/${sName}_fmdv_log.txt

vp1count=$((mapvp1))
cutoff=0
if [[ $vp1count -eq $cutoff ]]; then
	topvp1=""
fi
        
gencount=$((mapall))
if [[ $gencount -eq $cutoff ]]; then
    topgen="0-DEFAULT"
fi

stats=$(cat ${sPath}/${sName}_clip_depth_charge.txt | tail -n1 | sed 's/\t/,/g')


dep=$(echo $stats | cut -f3 -d ',')
echo "nanofmdv - average coverage = ${dep}" >> ${sPath}/${sName}_fmdv_log.txt

gen1=$(echo $stats | cut -f5 -d ',')
echo "nanofmdv - genome % covered with >= 1 reads = ${gen1}" >> ${sPath}/${sName}_fmdv_log.txt

gen5=$(echo $stats | cut -f7 -d ',')
echo "nanofmdv - genome % covered with >= 5 reads = ${gen5}" >> ${sPath}/${sName}_fmdv_log.txt

gen10=$(echo $stats | cut -f8 -d ',')
echo "nanofmdv - genome % covered with >= 10 reads = ${gen10}" >> ${sPath}/${sName}_fmdv_log.txt

gen20=$(echo $stats | cut -f9 -d ',')
echo "nanofmdv - genome % covered with >= 20 reads = ${gen20}" >> ${sPath}/${sName}_fmdv_log.txt

blast=$(cut -f2- ${sPath}/${sName}_fmdv_preconsensus_blast.txt | head -n1)
echo "nanofmdv - top BLAST hit of consensus = ${blast}" >> ${sPath}/${sName}_fmdv_log.txt

blast=$(head -n1 ${sPath}/${sName}_fmdv_preconsensus_blast.txt | cut -f2- | sed 's/\t/,/g')
blast_ref=$(head -n1 ${sPath}/${sName}_fmdv_preconsensus_blast.txt | cut -f2,2 | sed 's/,/_/g')

echo "${sName},${blast_ref},${reads},${mapall},${topgen},${mapvp1},${topvp1},${map},${stats},${blast}" >> ${sPath}/${sName}_fmdv_stats.csv

echo "nanofmdv - top BLAST hit = ${blast_ref}"
echo "nanofmdv - raw reads = ${reads}"

echo "nanofmdv - mapped reads = ${map} [${mapPer} %]"
echo "nanofmdv - average coverage = ${dep}"
echo "nanofmdv - % of ref genome covered by 1 or more reads = ${gen1}"
echo "nanofmdv - % of ref genome covered by 5 or more reads = ${gen5}"
echo "nanofmdv - % of ref genome covered by 10 or more reads = ${gen10}"
echo "nanofmdv - % of ref genome covered by 20 or more reads = ${gen20}"

if [ ! -n "$OPT_X" ]; then
	echo "nanofmdv - consensus length = ${conL} with non-N base length = ${conNonN}"
fi

#some cleanup options for after
#rm -rf ${sPath}/${sName}_vp1.bam ${sPath}/${sName}_vp1.bam.bai ${sPath}/${sName}_complete.bam ${sPath}/${sName}_complete.bam.bai

echo "nanofmdv - preconsensus length = ${preconL} with non-N base length = ${preconNonN}"  >> ${sPath}/${sName}_fmdv_log.txt
echo "nanofmdv - consensus length = ${conL} with non-N base length = ${conNonN}"  >> ${sPath}/${sName}_fmdv_log.txt