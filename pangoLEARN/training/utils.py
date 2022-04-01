from Bio import SeqIO
import os
import csv
import sys

import collections
import hashlib
import collections
import csv

def version_from_init(init_file):
    version=None
    with open(init_file, "r") as fr:
        for l in fr:
            if l.startswith("__version__"):
                l = l.rstrip("\n")
                version = l.split('=')[1]
                version = version.replace('"',"").replace(" ","")
                break
    return version

def get_pango_version(pango_path):
    version =""

    for r,d,f in os.walk(pango_path):
        for fn in f:
            if fn == "__init__.py":
                version = version_from_init(os.path.join(r, fn))
                if not version:
                    continue
    print("Pango version is:", version)

    if not version:
        sys.sterr.write("No version found at pango path")
        sys.exit(-1)
    else:
        return version

def get_hash_string(record):
    seq = str(record.seq).upper().encode()
    hash_object = hashlib.md5(seq)
    hash_string = hash_object.hexdigest()
    return hash_string

def get_dict(in_csv,name_column,data_column):
    this_dict = {}
    with open(in_csv,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            this_dict[row[name_column]] = row[data_column]
    return this_dict

def add_to_hash(seq_file):
    hash_map = {}
    seq_hash = {}
    for record in SeqIO.parse(seq_file, "fasta"):
        seq = str(record.seq).upper().encode()
        hash_object = hashlib.md5(seq)
        hash_map[hash_object.hexdigest()] = record.id
        seq_hash[str(record.seq)] = record.id
    return hash_map,seq_hash