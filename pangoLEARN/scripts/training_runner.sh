#!/bin/bash
source /home/s1680070/.bashrc_conda
conda activate grinch

TODAY=$(date +%F)
OUTDIR=data_release_$TODAY

mkdir /raid/shared/pangolearn_training/$OUTDIR

echo "Training model version: $TODAY"

LATEST_DATA=$(ls -td /raid/shared/2021* | head -n 1)

cd /raid/shared/pango-designation && git pull #gets any updates to the reports in the data directory
PANGO_V=$(git tag --sort=committerdate | tail -1)

echo "pango version $PANGO_V"

cd /raid/shared/pangolearn_training #gets any updates to the reports in the data directory

snakemake --snakefile /home/shared/pangoLEARN/pangoLEARN/scripts/curate_alignment.smk --nolock --configfile config.yaml --cores 1 --config outdir=$OUTDIR datadir=$LATEST_DATA pangolearn_version=$TODAY pango_version=$PANGO_V

cp /raid/shared/pangolearn_training/$OUTDIR/pangolearn.init.py /home/shared/pangoLEARN/pangoLEARN/__init__.py
cp /raid/shared/pangolearn_training/$OUTDIR/decision* /home/shared/pangoLEARN/pangoLEARN/data/
cp /raid/shared/pangolearn_training/$OUTDIR/metadata.downsample.csv /home/shared/pangoLEARN/pangoLEARN/data/lineages.downsample.csv
cp /raid/shared/pangolearn_training/$OUTDIR/lineage.hash.csv /home/shared/pangoLEARN/pangoLEARN/data/lineages.hash.csv

git add /home/shared/pangoLEARN/pangoLEARN/__init__.py
git add /home/shared/pangoLEARN/pangoLEARN/data/decision*
git add /home/shared/pangoLEARN/pangoLEARN/data/lineages*


