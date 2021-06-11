#!/usr/bin/env python3

import csv
import datetime
import operator
# import argparse
from Bio import SeqIO


# def parse_args():
#     parser = argparse.ArgumentParser(description="""Pick a representative sample for each unique sequence""",
#                                     formatter_class=argparse.RawTextHelpFormatter)
#     parser.add_argument('--in-metadata', dest = 'in_metadata', required=True, help='CSV of containing sequence_name and nucleotide_variants columns, the latter being | separated list of variants')
#     parser.add_argument('--in-fasta', dest = 'in_fasta', required=True, help='FASTA of all input seqs')
#     parser.add_argument('--diff', dest = 'diff', required=True, type=int, help='Samples within distance DIFF of included others may be excluded by the downsampler')
#     parser.add_argument('--out-metadata', dest = 'out_metadata', required=True, help='CSV to write out')
#     parser.add_argument('--out-fasta', dest = 'out_fasta', required=True, help='FASTA to write downsampled seqs')
#     parser.add_argument('--outgroups', dest = 'outgroups', required=False, help='Lineage splits file containing representative outgroups to protect')
#     parser.add_argument('--downsample_date_excluded', action='store_true', help='Downsample from those excluded as outside date window')
#     parser.add_argument('--downsample_included', action='store_true', help='Downsample from all included sequences')
#     parser.add_argument('--downsample_lineage_size', type=int, default=None, help='Min size of lineages to downsample, if unspecified no lineage-aware downsampling')

#     args = parser.parse_args()
#     return args

def parse_outgroups(outgroup_file):
    """
    input is CSV, last column being the representative outgroups:
    """
    outgroups = []
    if not outgroup_file:
        return outgroups
    with open(outgroup_file, "r") as outgroup_handle:
        line = outgroup_handle.readline()
        while line:
            try:
                outgroup = line.strip().split(",")[-1]
                outgroups.append(outgroup)
            except:
                continue
            line = outgroup_handle.readline()
    return(outgroups)

def get_count_dict(in_metadata):
    count_dict = {}
    num_samples = 0
    with open(in_metadata,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            num_samples += 1
            for var in row["nucleotide_variants"].split("|"):
                if var in count_dict:
                    count_dict[var] += 1
                else:
                    count_dict[var] = 1
    print("Found", len(count_dict), "variants")

    sorted_tuples = sorted(count_dict.items(), key=operator.itemgetter(1))
    count_dict = {k: v for k, v in sorted_tuples}
    return count_dict, num_samples

def get_lineage_dict(in_metadata, min_size):
    lineage_dict = {}
    if min_size is None:
        return lineage_dict

    with open(in_metadata,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if "lineage" in row:
                lin = row["lineage"]
                if lin in lineage_dict:
                    lineage_dict[lin].append(row["sequence_name"])
                else:
                    lineage_dict[lin] = [row["sequence_name"]]
    print("Found", len(lineage_dict), "lineages")

    small_lineages = [lin for lin in lineage_dict if len(lineage_dict[lin]) < min_size]
    for lin in small_lineages:
        del lineage_dict[lin]
    print("Found", len(lineage_dict), "lineages with at least", min_size, "representative sequences")

    return lineage_dict

def get_by_frequency(count_dict, num_samples, band=[0.1,1.0]):
    lower_bound = num_samples*band[0]
    upper_bound = num_samples*band[1]
    most_frequent = [k for k in count_dict if lower_bound < count_dict[k] <= upper_bound]
    print(len(most_frequent), "lie in frequency band", band)
    return most_frequent

def num_unique(muts1, muts2):
    u1 = [m for m in muts1 if m not in muts2]
    u2 = [m for m in muts2 if m not in muts1]
    return len(u1+u2)

def should_downsample_row(row, downsample_date_excluded=True, downsample_included=False, downsample_lineage_size=None, lineage_dict={}):
    if downsample_included and row["why_excluded"] in [None, "None", ""]:
        return True
    if downsample_date_excluded and row["why_excluded"] in [None, "None", ""] and "date_filter" in row and row["date_filter"].startswith("sample_date older than"):
        return True
    if downsample_lineage_size and row["lineage"] in lineage_dict:
        return True
    return False

def downsample(in_metadata, out_metadata, in_fasta, out_fasta, max_diff, outgroup_file, downsample_date_excluded, downsample_included, downsample_lineage_size):
    original_num_seqs = 0
    sample_dict = {}
    var_dict = {}

    count_dict, num_samples = get_count_dict(in_metadata)
    most_frequent = get_by_frequency(count_dict, num_samples, band=[0.05,1.0])
    very_most_frequent = get_by_frequency(count_dict, num_samples, band=[0.5,1.0])

    lineage_dict = get_lineage_dict(in_metadata,downsample_lineage_size)

    outgroups = parse_outgroups(outgroup_file)
    indexed_fasta = SeqIO.index(in_fasta, "fasta")

    with open(in_metadata, 'r', newline = '') as csv_in, \
        open(out_fasta, 'w', newline = '') as fa_out, \
        open(out_metadata, 'w', newline = '') as csv_out:

        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
        writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
        writer.writeheader()

        for row in reader:
            fasta_header = row["sequence_name"]
            if fasta_header not in indexed_fasta:
                continue
            if original_num_seqs % 1000 == 0:
                now = datetime.datetime.now()
                print("%s Downsampled from %i seqs to %i seqs" %(str(now), original_num_seqs, len(sample_dict)))
            original_num_seqs += 1

            if fasta_header in outgroups or not should_downsample_row(row,downsample_date_excluded, downsample_included,
                                                                      downsample_lineage_size,lineage_dict):
                if fasta_header in outgroups:
                    row["why_excluded"]=""
                writer.writerow(row)
                if row["why_excluded"] in [None, "None", ""] and fasta_header in indexed_fasta:
                    seq_rec = indexed_fasta[fasta_header]
                    fa_out.write(">" + seq_rec.id + "\n")
                    fa_out.write(str(seq_rec.seq) + "\n")
                else:
                    print(row["why_excluded"], fasta_header, (fasta_header in indexed_fasta))
                continue

            muts = row["nucleotide_variants"].split("|")
            if len(muts) < max_diff:
                #if not row["why_excluded"]:
                #    row["why_excluded"] = "downsampled with diff threshold %i" %max_diff
                writer.writerow(row)
                continue

            found_close_seq = False

            samples = set()
            low_frequency_muts = [mut for mut in muts if mut not in most_frequent]
            if len(low_frequency_muts) == 0:
                low_frequency_muts = [mut for mut in muts if mut not in very_most_frequent]
                if len(low_frequency_muts) == 0:
                    low_frequency_muts = muts
            if len(low_frequency_muts) > max_diff + 1:
                low_frequency_muts = low_frequency_muts[:max_diff+1]
            for mut in low_frequency_muts:
                if mut in var_dict:
                    samples.update(var_dict[mut])
            if downsample_lineage_size:
                samples = list( samples & set(lineage_dict[row["lineage"]]) )

            for sample in samples:
                if num_unique(muts, sample_dict[sample]) <= max_diff:
                    found_close_seq = True
                    #if not row["why_excluded"]:
                    #    row["why_excluded"] = "downsampled with diff threshold %i" %max_diff
                    writer.writerow(row)
                    break
            if not found_close_seq:
                sample_dict[fasta_header] = muts
                for mut in muts:
                    if mut not in var_dict:
                        var_dict[mut] = [fasta_header]
                    else:
                        var_dict[mut].append(fasta_header)
                row["why_excluded"] = ""
                writer.writerow(row)
                if fasta_header in indexed_fasta:
                    seq_rec = indexed_fasta[fasta_header]
                    fa_out.write(">" + seq_rec.id + "\n")
                    fa_out.write(str(seq_rec.seq) + "\n")

    now = datetime.datetime.now()
    print("%s Downsampled from %i seqs to %i seqs" %(str(now), original_num_seqs, len(sample_dict)))
    # return sample_dict.keys()

# def main():
#     args = parse_args()
#     subsample = downsample(args.in_metadata, args.out_metadata, args.in_fasta, args.out_fasta, args.diff, args.outgroups, args.downsample_date_excluded, args.downsample_included, args.downsample_lineage_size)

# if __name__ == '__main__':
#     main()