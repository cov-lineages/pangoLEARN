#!/bin/bash
source /localdisk/home/s1680070/.bashrc
conda activate pangolin



TODAY=$(date +%F)
OUTDIR=${TODAY}_pangoLEARN
echo $OUTDIR
REF="/localdisk/home/s1680070/repositories/pangolin/pangolin/data/reference.fasta"

mkdir /localdisk/home/shared/raccoon-dog/$OUTDIR

echo "Training model version: $TODAY"

LATEST_DATA=$(ls -td /localdisk/home/shared/raccoon-dog/2021*_gisaid | head -n 1)

cd /localdisk/home/s1680070/repositories/pango-designation && git pull #gets any updates to the reports in the data directory
PANGO_V=$(git tag --sort=committerdate | tail -1)
LINEAGES_CSV="/localdisk/home/s1680070/repositories/pango-designation/lineages.csv"

echo "pango version $PANGO_V"

cd /localdisk/home/shared/raccoon-dog/ #gets any updates to the reports in the data directory

snakemake --snakefile   /localdisk/home/s1680070/repositories/pangoLEARN/pangoLEARN/scripts/curate_alignment.smk --nolock --cores 1 --config lineages_csv=$LINEAGES_CSV reference=$REF =outdir=$OUTDIR datadir=$LATEST_DATA pangolearn_version=$TODAY pango_version=$PANGO_V

cp /localdisk/home/shared/raccoon-dog/$OUTDIR/pangolearn.init.py   /localdisk/home/s1680070/repositories/pangoLEARN/pangoLEARN/__init__.py
cp /localdisk/home/shared/raccoon-dog/$OUTDIR/decision*   /localdisk/home/s1680070/repositories/pangoLEARN/pangoLEARN/data/
cp /localdisk/home/shared/raccoon-dog/$OUTDIR/metadata.downsample.csv   /localdisk/home/s1680070/repositories/pangoLEARN/pangoLEARN/data/lineages.downsample.csv
cp /localdisk/home/shared/raccoon-dog/$OUTDIR/lineage.hash.csv   /localdisk/home/s1680070/repositories/pangoLEARN/pangoLEARN/data/lineages.hash.csv

git add   /localdisk/home/s1680070/repositories/pangoLEARN/pangoLEARN/__init__.py
git add   /localdisk/home/s1680070/repositories/pangoLEARN/pangoLEARN/data/decision*
git add   /localdisk/home/s1680070/repositories/pangoLEARN/pangoLEARN/data/lineages*


