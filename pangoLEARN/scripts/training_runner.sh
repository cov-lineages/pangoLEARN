#!/bin/bash
source /localdisk/home/s1680070/.bashrc
conda activate pangolin

TODAY=$(date +%F)
OUT=${TODAY}_pangoLEARN
OUTDIR=/localdisk/home/shared/raccoon-dog/$OUT

echo $OUTDIR

if [ -d $OUTDIR ] 
then
    echo "Directory $OUTDIR exists." 
else
    mkdir $OUTDIR 
    echo "Directory $OUTDIR does not exist, making it."
fi

echo "Training model version: $TODAY"

if [ -z "$1" ]
then
    DATA_DATE=$TODAY
else
    DATA_DATE=$1
fi

# LATEST_DATA=$(ls -td /localdisk/home/shared/raccoon-dog/2021*_gisaid/publish/gisaid | head -n 1)

cd /localdisk/home/s1680070/repositories/pango-designation && git pull #gets any updates to the reports in the data directory
PANGO_V=$(git tag --points-at HEAD)
echo "pango version $PANGO_V"

cd /localdisk/home/shared/raccoon-dog/ #gets any updates to the reports in the data directory
echo "--config outdir=$OUTDIR data_date=$DATA_DATE pangolearn_version=$TODAY pango_version=$PANGO_V"
echo "pangoLEARN training starting" | mail -s "update lineageTree.pb with pango designation version $PANGO_V" angie@soe.ucsc.edu
snakemake --snakefile /localdisk/home/s1680070/repositories/pangoLEARN/pangoLEARN/scripts/curate_alignment.smk --rerun-incomplete --nolock --cores 1 --config outdir=$OUTDIR data_date=$DATA_DATE pangolearn_version=$TODAY pango_version=$PANGO_V

cp $OUTDIR/pangolearn.init.py   /localdisk/home/s1680070/repositories/pangoLEARN/pangoLEARN/__init__.py
cp $OUTDIR/decision*   /localdisk/home/s1680070/repositories/pangoLEARN/pangoLEARN/data/
cp $OUTDIR/metadata.downsample.csv   /localdisk/home/s1680070/repositories/pangoLEARN/pangoLEARN/data/lineages.downsample.csv
cp $OUTDIR/lineage.hash.csv   /localdisk/home/s1680070/repositories/pangoLEARN/pangoLEARN/data/lineages.hash.csv

git add  /localdisk/home/s1680070/repositories/pangoLEARN/pangoLEARN/__init__.py
git add  /localdisk/home/s1680070/repositories/pangoLEARN/pangoLEARN/data/decision*
git add  /localdisk/home/s1680070/repositories/pangoLEARN/pangoLEARN/data/lineages*
