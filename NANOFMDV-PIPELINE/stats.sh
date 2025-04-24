#!/bin/bash

inDir=$1
echo "${inDir}"
nanoPath="${NANOFMDV_ENV}"

cat ${inDir}/*/*_fmdv_stats.csv > temp.csv

cat ${nanoPath}/stats_head.csv temp.csv > fmdv_all_stats.csv

rm -f temp.csv
