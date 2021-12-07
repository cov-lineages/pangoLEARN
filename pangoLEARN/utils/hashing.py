import os
import collections
import hashlib
import collections
import csv
from Bio import SeqIO
import pangoLEARN.utils.misc


def add_to_hash(seq_file):
    hash_map = {}
    seq_hash = {}
    for record in SeqIO.parse(seq_file, "fasta"):
        seq = str(record.seq).upper().encode()
        hash_object = hashlib.md5(seq)
        hash_map[hash_object.hexdigest()] = record.id
        seq_hash[str(record.seq)] = record.id
    return hash_map,seq_hash

def get_hash_string(record):
    seq = str(record.seq).upper().encode()
    hash_object = hashlib.md5(seq)
    hash_string = hash_object.hexdigest()
    return hash_string


def create_hash(designations,in_fasta,seq_hash,hashed_designations,hash_fasta):

    designated = misc.get_dict(designations,"sequence_name","lineage")
    hash_map,seq_hash_dict = add_to_hash(in_fasta)

    with open(out_hash,"w") as fw:
        fw.write("seq_hash,lineage\n")
        for seq_hash in hash_map:
            seq_name = hash_map[seq_hash]
            lin_name = designated[seq_name]
            fw.write(f"{seq_hash},{lin_name}\n")
    
    num_seqs = 0
    with open(hashed_designations, "w") as fw2:
        fw2.write("taxon,lineage\n")
        with open(hash_fasta, "w") as fw:
            for seq in seq_hash_dict:
                num_seqs +=1
                fw.write(f">{seq_hash_dict[seq]}\n{seq}\n")
                fw2.write(f"{seq_hash_dict[seq]},{designated[seq_hash_dict[seq]]}\n")

    print("Number of seqs going into training: ",f"{num_seqs}")