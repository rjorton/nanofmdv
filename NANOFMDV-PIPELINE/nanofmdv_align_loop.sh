#!/bin/bash


echo "NANOFMDV_ALIGN_LOOP SCRIPT STARTED"

while getopts "xn:r:t:p:m:d:xi:s:" TEST; do
	case $TEST in
	
	x) OPT_X=1
	;;
	n) nanoPath=$OPTARG
	;;
	r) refPath=$OPTARG
	;;
	t) threads=$OPTARG
	;;
	p) primersFile=$OPTARG
	;;
	m) medModel=$OPTARG
	;;
	d) depArg=$OPTARG
	;;

	i) inDir=$OPTARG
	;;
	s) samplesFile=$OPTARG
	;;
	
	esac
done



if [[ -n "${nanoPath}" ]]; then
	# remove trailing / if exists
	if [[ "$nanoPath" =~ '/'$ ]]; then 
		nanoPath=$(echo ${nanoPath} | sed 's/\/$//')
		
	fi
	
	echo "-n: NANOFMDV PATH SET TO = ${nanoPath}"

	if [ ! -d "$nanoPath" ]; then
		echo "ERROR - THE NANOFMDV DIRECTORY SUPPLIED VIA THE -n OPTION DOES NOT EXIST: ${nanoPath}"
		echo "EXITING"
		exit 1
    fi
else
	if [[ -n "${NANOFMDV_ENV}" ]]; then
		nanoPath="${NANOFMDV_ENV}"
		nanoPath=$(echo ${nanoPath} | sed 's/\/$//')

		echo "NANOFMDV_ENV PATH IS SET - NANOFMDV PATH = ${nanoPath}"

		if [ ! -d "$nanoPath" ]; then
			echo "ERROR - THE NANOFMDV DIRECTORY SUPPLIED VIA THE NANOFMDV_ENV DOES NOT EXIST: ${nanoPath}"
            echo "EXITING"
			exit 1
        fi
    else
		echo "ERROR - THERE IS NO NANOFMDV_ENV ENVIRONMENT VARIABLE SET OR PATH SUPPLIED"
		echo "SET A NANOFMDV ENV VARIABLE OR USE THE -n OPTION TO SPECIFY THE LOCATION OF THE NANOFMDV SCRIPTS FOLDER"
		echo "EXITING"
		exit 1
	fi
fi



if [[ -n "${refPath}" ]]; then
	refPath=$(echo ${refPath} | sed 's/\/$//')                
        
	echo "-r: NANOFMDV-REF PATH SET TO = ${refPath}"

	if [ ! -d "$refPath" ]; then
        echo "ERROR - THE NANOFMDV-REF DIRECTORY SUPPLIED VIA THE -r OPTION DOES NOT EXIST: ${refPath}"
		echo "EXITING"
		exit 1
	fi
else	
	refPath="${nanoPath}/Refs"
	echo "NANOFMDV-REF PATH SET TO = ${refPath}"

	if [ ! -d "$refPath" ]; then
        echo "ERROR - THE NANOFMDV-REF DIRECTORY (DERIVED FROM THE NANOFMDV PATH/ENV) DOES NOT EXIST: ${refPath}"
		echo "EXITING"
		exit 1
        fi
fi



if [ -n "$inDir" ]; then
	#if [[ "$inDir" =~ '/'$ ]]; then 
	#fi

	inDir=$(echo ${inDir} | sed 's/\/$//')
	echo "INPUT DIRECTORY = ${inDir}"
else
	inDir="."
	echo "NO INPUT DIRECTORY SET SO DEFAULTING TO: ${inDir}"
fi



if [ ! -d "$inDir" ]; then
	echo "ERROR - INPUT DIRECTOR DOES NOT EXIST: $inDir"
	echo "EXITING"
	exit 1
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






if [ -n "$samplesFile" ]
then
	if [ -f "$samplesFile" ]; then
        echo "CHECKING FORMAT OF SAMPLES-FILE: ${samplesFile} - INPUT DIR: ${inDir} - OUTPUT DIR: ${outDir}"
	else
		echo "ERROR - THE SAMPLES FILE SUPPLIED WITH THE -s OPTION DOES NOT EXIST: ${samplesFile}"
		echo "EXITING"
		exit 1
	fi

	sampleN=$(tail -n +2 ${samplesFile} | cut -f 2 -d ',' | wc -l | sed 's/ //g')
	sampleUN=$(tail -n +2 ${samplesFile} | cut -f 2 -d ',' | uniq | wc -l | sed 's/ //g')

	if [ "$sampleN" != "$sampleUN" ]; then
		echo "ERROR - DUPLICATE SAMPLE NAMES IN THE SAMPLES FILE - UNIQUE COUNTS BELOW:"
		echo "COUNT SAMPLE"
		tail -n +2 ${samplesFile} | cut -f 2 -d ',' | uniq -c | sort -n -k1,1
		echo "EXITING"
		exit 1
	fi

	barcodeN=$(tail -n +2 ${samplesFile} | cut -f 1 -d ',' | wc -l | sed 's/ //g')
	barcodeUN=$(tail -n +2 ${samplesFile} | cut -f 1 -d ',' | uniq | wc -l | sed 's/ //g')

	if [ "$barcodeN" != "$barcodeUN" ]; then
		echo "ERROR - DUPLICATE BARCODE NAMES IN THE SAMPLES FILE - UNIQUE COUNTS BELOW"
		echo "COUNT SAMPLE"
		tail -n +2 ${samplesFile} | cut -f 1 -d ',' | uniq -c | sort -n -k1,1
		echo "EXITING"
		exit 1
	fi

	sCount=0

	while read sLine; do
		#commaN=$(echo ${sLine} | sed 's/[^,]//g' | awk '{ print length }')
		commaN=$(echo ${sLine} | awk -F\, '{print NF-1}')
		if [ "$commaN" != "1" ]; then
			echo "ERROR - THERE ARE MORE THAN 2 FIELDS (1 COMMA) IN THE LINE: ${sLine}"
			echo "COMMA IS THE FIELD SEPARATOR - DO NOT USE A COMMA IN THE BARCODE OR SAMPLE NAMES"
			echo "EXITING"
			exit 1
		fi

		if [[ $sLine = *"/"* ]]; then
        	echo "ERROR - THERE ARE FORWARD SLASHES IN THE SAMPLES FILE LINE: ${sLine}"
            echo "EXITING"
            exit 1
        fi

		if [[ $sLine = *"\""* ]]; then
        	echo "ERROR - THERE ARE DOUBLE QUOTES IN THE SAMPLES FILE LINE: ${sLine}"
            echo "EXITING"
            exit 1
    	fi

		if [[ $sLine = *"'"* ]]; then
			echo "ERROR - THERE ARE SINGLE QUOTES IN THE SAMPLES FILE LINE: ${sLine}"
            echo "EXITING"
            exit 1
        fi

		#sBarcode=$(echo ${sLine} | cut -f1 -d ',' | xargs echo -n)
		sBarcode=$(echo ${sLine} | cut -f1 -d ',' | xargs)
		sName=$(echo ${sLine} | cut -f2 -d ',' | xargs)

		sCount=$((sCount + 1))

		if [[ "$sCount" -eq 1 ]]; then
			if [ "$sLine" != "barcode,sample" ]; then
 				echo "ERROR - HEADER LINE OF SAMPLES FILE IS NOT THE EXPECTED: barcode,sample"
				echo "EXITING"
				exit 1
			fi
		fi

		if [[ $sBarcode = *" "* ]]; then
            echo "ERROR - THERE ARE SPACES WITHIN THE BARCODE NAME: ${sBarcode}"
            echo "EXITING"
            exit 1
        fi

		if [[ $sName = *" "* ]]; then
            echo "ERROR - THERE ARE SPACES WITHIN THE SAMPLE NAME: ${sName}"
            echo "EXITING"
            exit 1
        fi

		if [[ "$sCount" -gt 1 ]]; then
			if [ ! -d "${inDir}/${sBarcode}" ]; then
				echo "ERROR - INPUT BARCODE SUBDIRECTORY DOES NOT EXIST: ${inDir}/${sBarcode}"
				echo "EXITING"
				exit 1

			fi
		fi

	done < "$samplesFile"

	sCount=$((sCount - 1))
	echo "SAMPLES IN SAMPLE FILE = ${sCount}"

	echo "ALIGNING EACH SAMPLE"
	sCount=0

	while read sLine; do

		sBarcode=$(echo ${sLine} | cut -f1 -d ',' | xargs)
		sName=$(echo ${sLine} | cut -f2 -d ',' | xargs)

		sCount=$((sCount + 1))

		if [ -n "$OPT_X" ]; then
			${nanoPath}/nanofmdv_single_barcode.sh -q ${inDir}/${sName}/${sName}.fastq -x -d ${depArg} -p ${primersFile} -n ${nanoPath} -r ${refPath}
		else
			${nanoPath}/nanofmdv_single_barcode.sh -q ${inDir}/${sName}/${sName}.fastq -d ${depArg} -p ${primersFile} -n ${nanoPath} -r ${refPath}
		fi


	done < "$samplesFile"


else
	echo "PROCESSING SAMPLES VIA ALL barcode* SUBDIRS IN THE INPUT DIR: ${inDir}"
	sCount=0

	for i in $(find ${inDir} -mindepth 1 -maxdepth 1 -type d | sort)
	do
        sPath=$(echo $i | sed 's/^.\///')
        sName=$(basename ${sPath})
        
        sCount=$((sCount + 1))
          
        if [ -n "$OPT_X" ]; then
			${nanoPath}/nanofmdv_single_barcode.sh -q ${inDir}/${sName}/${sName}.fastq -x -d ${depArg} -p ${primersFile} -n ${nanoPath} -r ${refPath}
		else
			${nanoPath}/nanofmdv_single_barcode.sh -q ${inDir}/${sName}/${sName}.fastq -d ${depArg} -p ${primersFile} -n ${nanoPath} -r ${refPath}
		fi

	done

	echo "BARCODE FOLDERS PROCESSED = ${sCount}"
fi

echo "NANOFMDV ALIGN LOOP SCRIPT FINISHED"



