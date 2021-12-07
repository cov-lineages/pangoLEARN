
def get_designated_seqs(designations_csv,
                        gisaid_metadata,
                        input_seqs,
                        matched_output,
                        not_matched,
                        matched_seqs):
    seqs_len = 0
    lineages = {}
    with open(designations_csv,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            lineages[row["taxon"]]=row["lineage"]
            
    matched = {}
    with open(matched_output,"w") as fw:
        with open(gisaid_metadata,"r") as f:

            reader = csv.DictReader(f)
            header = reader.fieldnames

            writer = csv.DictWriter(fw, fieldnames=header, lineterminator="\n")
            writer.writeheader()

            for row in reader:
                if name in lineages:
                    new_row = row
                    new_row["lineage"] = lineages[name]
                    writer.writerow(new_row)
                    matched[name]= 1
    
    not_matched_count = 0
    with open(not_matched,"w") as fw:
        fw.write("missed_taxon,lineage\n")
        for taxon in lineages:
            if taxon not in matched:
                not_matched_count +=1
                fw.write(f"{taxon},{lineage}\n")

    all_seqs = SeqIO.index(input_seqs,"fasta")

    written=0
    with open(matched_seqs,"w") as fw:
        for taxon in matched:
            try:
                record = all_seqs[taxon]
                fw.write(f">{taxon}\n{record.seq}\n")
                written +=1
            except:
                sys.stderr.write(f"Error: mismatch between gisaid metadata and sequence names. Taxon {taxon} not matched.\n")
                sys.exit(-1)

    print("Number of sequences going into pangolearn training",len(matched))
    print("Number of sequences found on gisaid written", written)

    print("Number of designations not matched", not_matched_count)

    if written != len(matched):
        sys.stderr.write(f"Error: mismatch between gisaid metadata and sequence names.\n")
        sys.exit(-1)