import os
import sys
sys.path.insert(0, '/localdisk/home/s1362711/staging/pangoLEARN')

from Bio import SeqIO

import collections
import hashlib
import collections
import csv
from pangoLEARN.training.get_lineage_positions import get_relevant_positions
from Bio import SeqIO
from datetime import date
today = date.today()

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

pangoLEARN_path = config["pangoLEARN_path"].rstrip("/")
pangolin_path = config["pangolin_path"].rstrip("/")
pango_designation_path = config["pango_designation_path"].rstrip("/")
quokka_path = config["quokka_path"].rstrip("/")

data_date = config["data_date"]
config["trim_start"] = 265
config["trim_end"] = 29674
config["lineages_csv"]=f"{pango_designation_path}/lineages.csv"
config["reference"] = f"{pangolin_path}/pangolin/data/reference.fasta"
config["outgroups"] = f"{pangoLEARN_path}/pangoLEARN/training/outgroups.csv"
config["genbank_ref"] = f"{pangoLEARN_path}/pangoLEARN/training/WH04.gb"
config["datadir"]= f"/localdisk/home/shared/raccoon-dog/{data_date}_gisaid/publish/gisaid"

rule all:
    input:
        os.path.join(config["outdir"],"alignment.filtered.fasta"),
        os.path.join(config["outdir"],"lineage_recall_report.txt"),
        os.path.join(config["outdir"],"pangolearn.init.py"),
        os.path.join(config["outdir"],"lineage.hash.csv")

rule make_init:
    output:
        init = os.path.join(config["outdir"],"pangolearn.init.py")
    run:
        pangolearn_new_v = config["pangolearn_version"]
        pango_version = config["pango_version"]
        with open(output.init,"w") as fw:
            fw.write(f'''_program = "pangoLEARN"
__version__ = "{pangolearn_new_v}"
PANGO_VERSION = "{pango_version}"
''')

rule filter_alignment:
    input:
        csv = config["lineages_csv"],
        fasta = os.path.join(config["datadir"],f"gisaid_{data_date}_all_alignment.fa"),
        full_csv = os.path.join(config["datadir"],f"gisaid_{data_date}_all_metadata.csv")
    output:
        fasta = os.path.join(config["outdir"],"alignment.filtered.fasta"),
        csv = os.path.join(config["outdir"],"lineages.metadata.filtered.csv"),
        csv_all = os.path.join(config["outdir"],"lineages.designated.csv")
    run:
        csv_len = 0
        seqs_len = 0
        lineages = {}
        all_lineages = {}
        all_len = 0
        with open(input.csv,"r") as f:
            for l in f:
                l = l.rstrip("\n")
                name,lineage = l.split(",")
                
                lineages[name]=lineage
                csv_len +=1
                all_lineages[name]=lineage
                all_len +=1
        with open(output.csv,"w") as fw:
            with open(output.csv_all,"w") as fw2:
                with open(input.full_csv,"r") as f:
                
                    reader = csv.DictReader(f)
                    header = reader.fieldnames

                    writer = csv.DictWriter(fw, fieldnames=header, lineterminator="\n")
                    writer.writeheader()

                    writer2 = csv.DictWriter(fw2, fieldnames=header, lineterminator="\n")
                    writer2.writeheader()

                    for row in reader:
                        name = row["sequence_name"].replace("SouthAfrica","South_Africa")
                        if name in lineages:
                            new_row = row
                            new_row["lineage"] = lineages[name]
                            writer.writerow(new_row)
                        if name in all_lineages:
                            all_row = row
                            all_row["lineage"] = all_lineages[name]
                            writer2.writerow(all_row)
        written = {}
        with open(output.fasta,"w") as fw:
            for record in SeqIO.parse(input.fasta, "fasta"):
                record.id = record.id.replace("SouthAfrica","South_Africa")
                if record.id in lineages and not record.id in written:
                    fw.write(f">{record.id}\n{record.seq}\n")

                    written[record.id]=1
                    seqs_len +=1
        
        print("Number of sequences in gisaid designated", all_len)
        print("Number of sequences going into pangolearn training",csv_len)
        print("Number of sequences found on gisaid", seqs_len)
        

rule align_with_minimap2:
    input:
        fasta = os.path.join(config["outdir"],"alignment.filtered.fasta"),
        reference = config["reference"]
    output:
        sam = os.path.join(config["outdir"],"alignment.sam")
    shell:
        """
        minimap2 -a -x asm5 -t {workflow.cores} \
        {input.reference:q} \
        {input.fasta:q} > {output.sam:q}
        """

rule get_variants:
    input:
        sam = os.path.join(config["outdir"],"alignment.sam")
    output:
        csv = os.path.join(config["outdir"],"variants.csv")
    shell:
        """
        gofasta sam variants -t {workflow.cores} \
        --samfile {input.sam:q} \
        --reference {config[reference]} \
        --genbank {config[genbank_ref]} \
        --outfile {output.csv}
        """

rule add_lineage:
    input:
        csv = os.path.join(config["outdir"],"variants.csv"),
        lineages = os.path.join(config["outdir"],"lineages.designated.csv")
    output:
        csv = os.path.join(config["outdir"],"variants.lineages.csv")
    run:
        lineages_dict = {}
        with open(input.lineages,"r") as f:
            reader= csv.DictReader(f)
            for row in reader:
                lineages_dict[row["sequence_name"]] = row["lineage"]
        with open(output.csv, "w") as fw:
            fw.write("sequence_name,nucleotide_variants,lineage,why_excluded\n")
            with open(input.csv,"r") as f:
                for l in f:
                    l = l.strip("\n")
                    name,variants = l.split(",")
                    if name =="query":
                        pass
                    elif name in lineages_dict:
                        fw.write(f"{name},{variants},{lineages_dict[name]},\n")
                        
rule filter_metadata:
    input:
        csv = os.path.join(config["outdir"],"variants.lineages.csv"),
        fasta = os.path.join(config["outdir"],"alignment.filtered.fasta")
    output:
        csv = os.path.join(config["outdir"],"metadata.final.csv")
    run:
        in_list = {}
        for record in SeqIO.parse(input.fasta,"fasta"):
            in_list[record.id] = 1

        with open(output.csv, "w") as fw:
            fw.write("sequence_name,lineage\n")
            with open(input.csv, "r") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    if row["sequence_name"] in in_list:
                        name = row["sequence_name"]
                        lineage = row["lineage"]
                        fw.write(f"{name},{lineage}\n")

rule get_relevant_postions:
    input:
        fasta = os.path.join(config["outdir"],"alignment.filtered.fasta"),
        csv = os.path.join(config["outdir"],"metadata.final.csv"),
        reference = config["reference"]
    output:
        relevant_pos_obj = os.path.join(config["outdir"],"relevantPositions.pickle"),
    run:
        get_relevant_positions(input.csv,input.fasta,input.reference,output.relevant_pos_obj)

rule run_training:
    input:
        fasta = os.path.join(config["outdir"],"alignment.filtered.fasta"),
        csv = os.path.join(config["outdir"],"metadata.final.csv"),
        reference = config["reference"],
        relevant_pos_obj = rules.get_relevant_postions.output.relevant_pos_obj
    params:
        path_to_script = pangoLEARN_path
    output:
        headers = os.path.join(config["outdir"],"randomForestHeaders_v1.joblib"),
        model = os.path.join(config["outdir"],"randomForest_v1.joblib"),
        txt = os.path.join(config["outdir"],"training_summary.txt")
    shell:
        """
        python {params.path_to_script}/pangoLEARN/training/pangoLEARNRandomForest_v1.py \
        {input.csv:q} \
        {input.fasta} \
        {input.reference:q} \
        {config[outdir]} \
        {input.relevant_pos_obj} \
        > {output.txt:q}
        """

rule get_recall:
    input:
        txt = rules.run_training.output.txt
    params:
        path_to_script = pangoLEARN_path
    output:
        txt = os.path.join(config["outdir"],"lineage_recall_report.txt")
    shell:
        """
        python {params.path_to_script}/pangoLEARN/training/processOutputFile.py {input.txt} > {output.txt}
        """

rule create_hash:
    input:
        fasta = os.path.join(config["outdir"],"alignment.filtered.fasta"),
        lin_designation = os.path.join(config["outdir"],"lineages.designated.csv")
    output:
        csv = os.path.join(config["outdir"],"lineage.hash.csv"),
        fasta = os.path.join(config["outdir"],"lineage.hash.fasta"),
        hashed_designations = os.path.join(config["outdir"],"designations.hash.csv")
    run:
        designated = get_dict(input.lin_designation,"sequence_name","lineage")

        hash_map,seq_hash_dict = add_to_hash(input.fasta)

        with open(output.csv,"w") as fw:
            fw.write("seq_hash,lineage\n")
            for seq_hash in hash_map:
                seq_name = hash_map[seq_hash]
                set_name = designated[seq_name]
                fw.write(f"{seq_hash},{set_name}\n")
        
        num_seqs = 0
        with open(output.hashed_designations, "w") as fw2:
            fw2.write("taxon,lineage\n")
            with open(output.fasta, "w") as fw:
                for seq in seq_hash_dict:
                    num_seqs +=1
                    fw.write(f">{seq_hash_dict[seq]}\n{seq}\n")
                    fw2.write(f"{seq_hash_dict[seq]},{designated[seq_hash_dict[seq]]}\n")

        print("Number of seqs going into training: ",f"{num_seqs}")
