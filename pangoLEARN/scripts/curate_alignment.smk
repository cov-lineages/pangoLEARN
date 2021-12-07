
import csv
from Bio import SeqIO
import os
import collections
import hashlib
import collections
import csv
from Bio import SeqIO
from pangoLEARN.training import downsample
from pangoLEARN.training.get_lineage_positions import get_relevant_positions
from pangoLEARN.training import process_designations

import pangoLEARN.utils.config as cfg
import pangoLEARN.utils.hashing
import pangoLEARN.utils.misc

from datetime import date
today = date.today()

cfg.min_config(config)

cfg.setup_paths(config)

rule all:
    input:
        os.path.join(config["outdir"],"alignment.filtered.fasta"),
        os.path.join(config["outdir"],"decision_tree_rules.zip"),
        os.path.join(config["outdir"],"pangolearn.init.py"),
        os.path.join(config["outdir"],"lineage.hash.csv")

rule make_init:
    output:
        init = os.path.join(config["outdir"],"pangolearn.init.py")
    run:
        create_init(config["pangoLEARN_version"],config["pango_version"],output.init)

rule filter_alignment:
    input:
        csv = config["lineages_csv"],
        fasta = os.path.join(config["datadir"],f"gisaid_{data_date}_all_alignment.fa"),
        full_csv = os.path.join(config["datadir"],f"gisaid_{data_date}_all_metadata.csv")
    output:
        matched_seqs = os.path.join(config["outdir"],"alignment.filtered.fasta"),
        not_matched = os.path.join(config["outdir"],"designated.not_matched.csv"),
        lineages_csv = os.path.join(config["outdir"],"lineages.designated.csv")
    run:
        process_designations.get_designated_seqs(input.csv,
                        input.full_csv,
                        input.fasta,
                        output.lineages_csv,
                        output.not_matched,
                        output.matched_seqs)
        
rule align_with_minimap2:
    input:
        fasta = rules.filter_alignment.output.matched_seqs,
        reference = config["reference"]
    output:
        sam = os.path.join(config["outdir"],"alignment.sam")
    shell:
        """
        minimap2 -a -x asm20 -t {workflow.cores} \
        {input.reference:q} \
        {input.fasta:q} > {output.sam:q}
        """

rule get_variants:
    input:
        sam = rules.align_with_minimap2.output.sam,
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

rule add_lineage_to_variants:
    input:
        csv = rules.get_variants.output.csv,
        lineages = rules.filter_alignment.output.lineages_csv
    output:
        csv = os.path.join(config["outdir"],"variants.lineages.csv")
    run:
        lineages_dict = misc.get_dict(input.lineages,"sequence_name","lineage")

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
        csv = rules.add_lineage_to_variants.output.csv,
        fasta = rules.filter_alignment.output.matched_seqs
    output:
        csv = os.path.join(config["outdir"],"metadata.copy.csv"),
        fasta = os.path.join(config["outdir"],"alignment.downsample.fasta")
    run:
        downsample.downsample(
            input.csv, 
            output.csv, 
            input.fasta, 
            output.fasta, 
            1, config["outgroups"], 
            False, 
            False, 
            10)

rule filter_metadata_to_just_downsample:
    input:
        csv = rules.filter_alignment.output.lineages_csv,
        fasta = rules.downsample.output.fasta
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

rule get_relevant_postions:
    input:
        fasta = rules.downsample.output.fasta,
        csv = rules.filter_metadata_to_just_downsample.output.csv,
        reference = config["reference"]
    output:
        relevant_pos_obj = os.path.join(config["outdir"],"relevantPositions.pickle"),
    run:
        get_relevant_positions(input.csv,input.fasta,input.reference,output.relevant_pos_obj)

rule run_training:
    input:
        fasta = rules.downsample.output.fasta,
        csv = rules.filter_metadata_to_just_downsample.output.csv,
        reference = config["reference"],
        relevant_pos_obj = rules.get_relevant_postions.output.relevant_pos_obj
    params:
        path_to_script = pangoLEARN_path
    output:
        headers = os.path.join(config["outdir"],"decisionTreeHeaders_v1.joblib"),
        model = os.path.join(config["outdir"],"decisionTree_v1.joblib"),
        txt = os.path.join(config["outdir"],"training_summary.txt")
    shell:
        """
        python {params.path_to_script}/pangoLEARN/training/pangoLEARNDecisionTree_v1.py \
        {input.csv:q} \
        {input.fasta} \
        {input.reference:q} \
        {config[outdir]} \
        {input.relevant_pos_obj} \
        A \
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

rule get_decisions:
    input:
        headers = os.path.join(config["outdir"],"decisionTreeHeaders_v1.joblib"),
        model = os.path.join(config["outdir"],"decisionTree_v1.joblib"),
        txt = rules.run_training.output.txt
    params:
        path_to_script = pangoLEARN_path
    output:
        txt = os.path.join(config["outdir"],"tree_rules.txt"),
        zipped = os.path.join(config["outdir"],"decision_tree_rules.zip")
    shell:
        """
        python {params.path_to_script}/pangoLEARN/training/getDecisionTreeRules.py \
        {input.model:q} {input.headers:q} {input.txt:q} \
        > {output.txt:q} && zip {output.zipped:q} {output.txt:q}
        """

rule create_hash:
    input:
        fasta = os.path.join(config["outdir"],"alignment.filtered.fasta"),
        lin_designation = os.path.join(config["outdir"],"lineages.designated.csv")
    output:
        seq_hash = os.path.join(config["outdir"],"lineage.hash.csv"),
        fasta = os.path.join(config["outdir"],"lineage.hash.fasta"),
        hashed_designations = os.path.join(config["outdir"],"designations.hash.csv")
    run:
        hashing.create_hash(input.lineage_designations,input.fasta,output.seq_hash,output.hashed_designations,output.fasta)
