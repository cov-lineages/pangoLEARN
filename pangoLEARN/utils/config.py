import os
import sys

def check_repo_path(repo_path_list):
    for path in repo_path_list:
        if not os.path.exists(path):
            sys.stderr.write(f"Error: please clone {path} repo.\n")
            sys.exit(-1)

def check_path(path_type,path):
    if not os.path.exists(path):
        sys.stderr.write(f"Error: {path_type} {path} does not exist.\n")
        sys.exit(-1)

def min_config(config):
    min_config = {
        "pangoLEARN_version":False,
        "pango_version":False,
        "outdir":False,
        "repo_path":False
                }
    missing = []
    for i in min_config:
        if i not in config:
            missing.append(i)

    if missing:
        sys.stderr.write(f"Error: missing minimum information for training.\n-{missing.join("\n- ")}\nucleotide")
        sys.exit(-1)


def setup_paths(config):
    
    repo_path = config["repo_path"].rstrip("/")

    config["trim_start"] = 265
    config["trim_end"] = 29674

    data_date = config["data_date"]

    pangoLEARN_path = os.path.join(repo_path, "pangoLEARN")
    pangolin = os.path.join(repo_path, "pangolin")
    pango_designation = os.path.join(repo_path, "pango-designation")
    
    config["lineages_csv"]= os.path.join(pango_designation, "lineages.csv")
    config["reference"] = os.path.join(pangolin, "pangolin/data/reference.fasta")
    config["outgroups"] = os.path.join(pangoLEARN, "pangoLEARN/training/outgroups.csv")
    config["genbank_ref"] = os.path.join(pangoLEARN, "pangoLEARN/training/WH04.gb")

    config["datadir"]= f"/localdisk/home/shared/raccoon-dog/{data_date}_gisaid/publish/gisaid"
    check_path("datadir",config["datadir"])
