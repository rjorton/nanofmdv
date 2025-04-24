#!/bin/bash


# single sample - NANOFMDV env or path - input fastq - use as stub
# loop samples - over all sample

# this script simply concatenates multiple FASTQs into single FASTQs
# either by looping over all barcode folders in a given input folder (default is current folder .)
# or by supplying a sample.csv file (barcode,sample) specifying barcode and sample names
# output is to NANOFMDV-ANALYSIS-YYYYMMDD or override with -o
# i = inputDir
# o = outputDir
# g = mac gzcat instead of zcat
# u = unzipped FASTQs
# s = samples.csv file
# n = location of NANOFMDV scripts folder - will assign Refs to Refs subfolder
# r = location of NANOFMDV Refs folder - overrides -p
# set a NANOFMDV_ENV OPTION by:
# export NANOFMDV_ENV=/path/to/nanofmdv/Scripts

# ToDo
# REMEMBER TO UPDATE PRIMER GOPRIME SCRIPT
#



echo "NANOFMDV PREPARE_READS SCRIPT STARTED"

while getopts "guvi:o:s:n:r:" TEST; do
	case $TEST in

    g) OPT_G=1
    ;;
	u) OPT_U=1
    ;;
	v) OPT_V=1
	;;
    i) inDir=$OPTARG
	;;
	o) outDir=$OPTARG
	;;
	s) samplesFile=$OPTARG
	;;
	n) nanoPath=$OPTARG
	;;
	r) refPath=$OPTARG
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



outDate=$(date +%Y%m%d)
outStub="NANOFMDV-ANALYSIS"

if [ -n "$outDir" ]; then
	outDir=$(echo ${outDir} | sed 's/\/$//')
    echo "OUTPUT DIRECTORY = ${outDir}"
else
    outDir="${outStub}-${outDate}"
	echo "NO OUTPUT DIRECTORY SET SO DEFAULTING TO: ${outDir}"
fi

if [ -d "$outDir" ]; then
	if [ -n "$OPT_V" ]; then
		echo "-v OPTION - OUTPUT DIRECT ALREADY EXISTS BUT OVERWRITING: ${outDir}"
	else
		echo "ERROR - EXITING AS OUTPUT DIRECTOR ALREADY EXISTS: $outDir"
		echo "TO OVERWRITE EXISTING RESULTS USE THE OPTION:  -v"
		echo "EXITING"
		exit 1
	fi
else
        mkdir ${outDir}

	if [ -d "$outDir" ]; then
		echo "CREATED OUTPUT DIRECTORY: $outDir"
	else
		echo "ERROR - COULD NOT CREATE THE OUTPUT DIRECTORY: $outDir"
		echo "CHECK PATHS AND PERMISSIONS"
		exit 1
	fi
fi



if [ -n "$OPT_U" ]; then
	echo "-u OPTION IS SET SO SETTING INPUT TO UNZIPPED FASTQ FILES ENDING WITH .fastq OR .fq"
fi



if [ -n "$OPT_G" ]; then
    echo "-g OPTION IS SET SO CHANGING ZCAT TO GZCAT FOR DECOMPRESSION OF FASTQ FILES ENDING WITH .fastq.gz or .fq.gz"
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



	echo "PROCESSING READS FOR EACH SAMPLE"
	sCount=0

	while read sLine; do

		sBarcode=$(echo ${sLine} | cut -f1 -d ',' | xargs)
		sName=$(echo ${sLine} | cut -f2 -d ',' | xargs)

		sCount=$((sCount + 1))

		if [[ "$sCount" -gt 1 ]]; then
			echo "BARCODE SUBFOLDER = ${sBarcode}, OUTPUT SUBFOLDER = ${sName}"

			mkdir -p ${outDir}/${sName}

            if [ -n "$OPT_U" ]; then
            	cat ${inDir}/${sBarcode}/*.f*q > ${outDir}/${sName}/${sName}.fastq
            else
                if [ -n "$OPT_G" ]; then
                    gzcat ${inDir}/${sBarcode}/*.f*q.gz > ${outDir}/${sName}/${sName}.fastq
                else
                    zcat ${inDir}/${sBarcode}/*.f*q.gz > ${outDir}/${sName}/${sName}.fastq
                fi
            fi
		fi

	done < "$samplesFile"


else
	echo "PROCESSING SAMPLES VIA ALL barcode* SUBDIRS IN THE INPUT DIR: ${inDir}"
	sCount=0

	for i in $(find ${inDir} -mindepth 1 -maxdepth 1 -type d)
	do
        sPath=$(echo $i | sed 's/^.\///')
        sName=$(basename ${sPath})
          
        if [[ $sName == barcode* ]]; then
            echo "BARCODE SUBFOLDER = ${sName}"
			sCount=$((sCount + 1))

            mkdir -p ${outDir}/${sName}

            if [ -n "$OPT_U" ]
            then
                cat ${inDir}/${sName}/*.f*q > ${outDir}/${sName}/${sName}.fastq
            else
                if [ -n "$OPT_G" ]; then
                	gzcat ${inDir}/${sName}/*.f*q.gz > ${outDir}/${sName}/${sName}.fastq
                else
                    zcat ${inDir}/${sName}/*.f*q.gz > ${outDir}/${sName}/${sName}.fastq
                fi
            fi
        fi
	done

	echo "BARCODE FOLDERS PROCESSED = ${sCount}"
fi

echo "NANOFMDV PREPARE READS SCRIPT FINISHED"
