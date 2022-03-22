#!/bin/bash

#source /localdisk/home/s1680070/.bashrc
#conda activate pangolin

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
    UID="s1362711"
else
    UID=$1
fi

if [ -z "$2" ]
then
    DATA_DATE=$TODAY
else
    DATA_DATE=$2
fi
echo "$UID training version $DATA_DATE"

# LATEST_DATA=$(ls -td /localdisk/home/shared/raccoon-dog/2021*_gisaid/publish/gisaid | head -n 1)

REPO_PATH=/localdisk/home/s1362711/staging

PANGO_PATH=$REPO_PATH/pango-designation
PLEARN_PATH=$REPO_PATH/pangoLEARN
PANGOLIN_PATH=$REPO_PATH/pangolin
QUOKKA_PATH=$REPO_PATH/quokka

echo "pango designation path $PANGO_PATH"
echo "pangoLEARN path $PLEARN_PATH"
echo "pangolin path $PANGOLIN_PATH"


cd $PANGO_PATH && git pull #gets any updates to the reports in the data directory
PANGO_V=$(git tag --points-at HEAD)
echo "pango version $PANGO_V"

cd /localdisk/home/shared/raccoon-dog/ #gets any updates to the reports in the data directory
echo "--config outdir=$OUTDIR data_date=$DATA_DATE pangolearn_version=$TODAY pango_version=$PANGO_V"
echo "pangoLEARN training starting" | mail -s "update lineageTree.pb with pango designation version $PANGO_V" angie@soe.ucsc.edu
snakemake --snakefile $PLEARN_PATH/pangoLEARN/scripts/curate_alignment.smk --rerun-incomplete --nolock --cores 1 --config pango_designation_path=$PANGO_PATH pangolin_path=$PANGOLIN_PATH pangoLEARN_path=$PLEARN_PATH quokka_path=$QUOKKA_PATH outdir=$OUTDIR data_date=$DATA_DATE pangolearn_version=$TODAY pango_version=$PANGO_V

# cp $OUTDIR/pangolearn.init.py   /localdisk/home/s1680070/repositories/pangoLEARN/pangoLEARN/__init__.py
# cp $OUTDIR/decision*   /localdisk/home/s1680070/repositories/pangoLEARN/pangoLEARN/data/
# cp $OUTDIR/metadata.downsample.csv   /localdisk/home/s1680070/repositories/pangoLEARN/pangoLEARN/data/lineages.downsample.csv
# cp $OUTDIR/lineage.hash.csv   /localdisk/home/s1680070/repositories/pangoLEARN/pangoLEARN/data/lineages.hash.csv
