#!/bin/bash

inDir=$1
echo "${inDir}"
nanoPath="${NANOFMDV_ENV}"

cat ${inDir}/*/*_fmdv_stats.csv | grep -v "^Sample"> temp.csv

cat ${nanoPath}/stats_head.csv temp.csv | cut -f 1-3,6- -d ',' > fmdv_all_stats.csv

rm -f temp.csv


