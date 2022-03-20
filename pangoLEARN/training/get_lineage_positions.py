import csv
import pickle
import collections
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align.AlignInfo import SummaryInfo

## outputs a dataframe of positions to input in training

def get_lineage_cns50_sites(lineage, lineage_taxa, sequences_index, reference_sequence): 
            
    ### output of list of sequence_IDs for the given lineage
    
    lineage_seqs = MultipleSeqAlignment([])
    for taxon in lineage_taxa:
        lineage_seqs.append(sequences_index[taxon])

    info = SummaryInfo(lineage_seqs)
    consensus_sequence =  info.gap_consensus(
    threshold=0.50, 
    ambiguous='N')
    
    nuc_position = []
    for i in range(len(consensus_sequence)):
        if reference_sequence[i] != "N":
            if consensus_sequence[i] != "N":
                if consensus_sequence[i] != reference_sequence[i]:
                    nuc_position.append(i)
    return nuc_position

def get_relevant_positions(designation_file,seq_file,ref_file,outfile):
    reference = SeqIO.read(ref_file, "fasta")
    sequences_index = SeqIO.index(seq_file, "fasta")
    lineage_designations = collections.defaultdict(list)
    lineage_set = set()

    with open(designation_file, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            lineage_designations[row["lineage"]].append(row["sequence_name"])
            lineage_set.add(row["lineage"])
    
    final_positions = set()
    for lineage in lineage_set:
        print(f"Getting positions for lineage {lineage}")
        positions = get_lineage_cns50_sites(lineage, lineage_designations[lineage], sequences_index, reference)
        print(f"\tFound {len(positions)}")
        for i in positions:
            final_positions.add(i)

    with open(outfile, 'wb') as pickle_file:
        pickle.dump(final_positions, pickle_file)