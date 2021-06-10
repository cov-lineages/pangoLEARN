import csv
from Bio import SeqIO
import os
import collections


config["trim_start"] = 265
config["trim_end"] = 29674
config["lineages_of_concern"] = ["B.1.1.7","B.1.351","P.1","P.2"]
config["outgroups"] = "/raid/shared/grinch/grinch/data/outgroups.csv"
config["genbank_ref"] = "/raid/shared/grinch/grinch/data/WH04.gb"

rule all:
    input:
        os.path.join(config["outdir"],"alignment.filtered.fasta"),
        os.path.join(config["outdir"],"decision_tree_rules.txt"),
        os.path.join(config["outdir"],"pangolearn.init.py")

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
        fasta = os.path.join(config["datadir"],"0","gisaid.UH.RD.aligned.fasta"),
        full_csv = os.path.join(config["datadir"],"2","lineages.metadata.csv")
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
                        

rule downsample:
    input:
        csv = os.path.join(config["outdir"],"variants.lineages.csv"),
        fasta = os.path.join(config["outdir"],"alignment.filtered.fasta")
    output:
        csv = os.path.join(config["outdir"],"metadata.copy.csv"),
        fasta = os.path.join(config["outdir"],"alignment.downsample.fasta")
    shell:
        """
        python /raid/shared/grinch/grinch/scripts/downsample.py --in-fasta {input.fasta} \
        --in-metadata {input.csv} \
        --downsample_lineage_size 10 \
        --diff 1 \
        --out-fasta {output.fasta} \
        --outgroups {config[outgroups]} \
        --out-metadata {output.csv}
        """

rule filter_metadata:
    input:
        csv = os.path.join(config["outdir"],"metadata.copy.csv"),
        fasta = os.path.join(config["outdir"],"alignment.downsample.fasta")
    output:
        csv = os.path.join(config["outdir"],"metadata.downsample.csv")
    run:
        in_downsample = {}
        for record in SeqIO.parse(input.fasta,"fasta"):
            in_downsample[record.id] = 1

        with open(output.csv, "w") as fw:
            fw.write("sequence_name,lineage\n")
            with open(input.csv, "r") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    if row["sequence_name"] in in_downsample:
                        name = row["sequence_name"]
                        lineage = row["lineage"]
                        fw.write(f"{name},{lineage}\n")

rule run_training:
    input:
        fasta = os.path.join(config["outdir"],"alignment.downsample.fasta"),
        csv = os.path.join(config["outdir"],"metadata.downsample.csv"),
        reference = config["reference"]
    output:
        headers = os.path.join(config["outdir"],"decisionTreeHeaders_v1.joblib"),
        model = os.path.join(config["outdir"],"decisionTree_v1.joblib"),
        txt = os.path.join(config["outdir"],"training_summary.txt")
    shell:
        """
        python /raid/shared/grinch/grinch/scripts/pangoLEARNDecisionTree_v1.py \
        {input.csv:q} \
        {input.fasta} \
        {input.reference:q} \
        {config[outdir]} \
        > {output.txt:q}
        """

rule get_recall:
    input:
        txt = rules.run_training.output.txt
    output:
        txt = os.path.join(config["outdir"],"lineage_recall_report.txt")
    shell:
        """
        python /raid/shared/grinch/grinch/scripts/processOutputFile.py {input.txt} > {output.txt}
        """

rule get_decisions:
    input:
        headers = os.path.join(config["outdir"],"decisionTreeHeaders_v1.joblib"),
        model = os.path.join(config["outdir"],"decisionTree_v1.joblib"),
        txt = rules.run_training.output.txt
    output:
        txt = os.path.join(config["outdir"],"decision_tree_rules.txt")
    shell:
        """
        python /raid/shared/grinch/grinch/scripts/getDecisionTreeRules.py \
        {input.model:q} {input.headers:q} {input.txt:q} \
        > {output.txt:q}
        """
