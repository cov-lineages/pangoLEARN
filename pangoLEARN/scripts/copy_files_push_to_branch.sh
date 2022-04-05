#!/bin/bash

OUTDIR=$1
PANGO_VERSION=$2
PLEARN_VERSION=$3
REPO_PATH=/localdisk/home/s1680070/repositories

cd $REPO_PATH/pangoLEARN && git pull

cp $OUTDIR/pangolearn.init.py   $REPO_PATH/pangoLEARN/pangoLEARN/__init__.py
cp $OUTDIR/pangolin_data.init.py   $REPO_PATH/pangolin-data/pangolin_data/__init__.py

cp $OUTDIR/decisionTreeHeaders_v1.joblib   $REPO_PATH/pangoLEARN/pangoLEARN/data/decisionTreeHeaders_v1.joblib
cp $OUTDIR/decisionTree_v1.joblib   $REPO_PATH/pangoLEARN/pangoLEARN/data/decisionTree_v1.joblib
cp $OUTDIR/decision_tree_rules.zip   $REPO_PATH/pangoLEARN/pangoLEARN/data/decision_tree_rules.zip

cp $OUTDIR/random*   $REPO_PATH/pangoLEARN/pangoLEARN/data/
cp $OUTDIR/random*   $REPO_PATH/pangoLEARN/pangoLEARN/data/

cp $OUTDIR/metadata.final.csv   $REPO_PATH/pangoLEARN/pangoLEARN/data/lineages.downsample.csv
cp $OUTDIR/lineage.hash.csv   $REPO_PATH/pangoLEARN/pangoLEARN/data/lineages.hash.csv

cp $OUTDIR/random*   $REPO_PATH/pangolin-data/pangolin_data/data/
cp $OUTDIR/random*   $REPO_PATH/pangolin-data/pangolin_data/data/
cp $OUTDIR/lineage.hash.csv   $REPO_PATH/pangolin-data/pangolin_data/data/lineages.hash.csv
cp $REPO_PATH/pango-designation/pango_designation/alias_key.json   $REPO_PATH/pangolin-data/pangolin_data/data/

echo $PLEARN_VERSION
cd $REPO_PATH/pangoLEARN && git checkout -b "prerelease_$PLEARN_VERSION"
git add $REPO_PATH/pangoLEARN/pangoLEARN/data/*
git add $REPO_PATH/pangoLEARN/pangoLEARN/__init__.py
git commit -m "adding latest decision tree and rf model to pangoLEARN repo for trained version $PLEARN_VERSION corresponding to $PANGO_VERSION"
git push --set-upstream origin "prerelease_$PLEARN_VERSION"
git checkout master

cd $REPO_PATH/pangolin-data && git checkout -b "prerelease_$PANGO_VERSION"
git add $REPO_PATH/pangolin-data/pangolin_data/data/*
git add $REPO_PATH/pangolin-data/pangolin_data/__init__.py
git commit -m "adding latest hash, alias file and rf model corresponding to $PANGO_VERSION"
git push --set-upstream origin "prerelease_$PANGO_VERSION"
git checkout main
