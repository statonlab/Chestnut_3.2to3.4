#!/data/apps/python/3.2.1/bin/python
import re, sys, getopt
from Bio import SeqIO

# This script will rename sequences from the updated annotation to match the original, and then replace the old predicted genes with the new ones.
# Note that for genes that were split, you will find that the input tsv file contains a third column with the original gene name, and a column with an altered original gene name.

input_tsv = "/staton/projects/chestnut/psudochro/analysis_081718_annotation/12k_RNA_Qrobur_081718/8_fixInternalStops/3_renameGenes/split_combined_genes.txt"

old_annotations = "/staton/projects/chestnut/psudochro/analysis_081718_annotation/12k_RNA_Qrobur_081718/8_fixInternalStops/3_renameGenes/Castanea_mollissima_scaffolds_v3.2_HQcds_iter1.fna"

new_annotations = "/staton/projects/chestnut/psudochro/analysis_081718_annotation/12k_RNA_Qrobur_081718/8_fixInternalStops/3_renameGenes/augustus.hints.codingseq"

final_annotations = "/staton/projects/chestnut/psudochro/analysis_081718_annotation/12k_RNA_Qrobur_081718/8_fixInternalStops/3_renameGenes/Castanea_mollissima_scaffolds_v3.4_HQcds_iter2.fna"

matches = {}
updated_names = {}
with open(input_tsv) as t:
    for line in t:
        line = line.rstrip().split()
        if len(line) == 2:
            oldName = line[0]
            newName = line[1]
            matches[newName] = oldName
        elif len(line) == 3:
            oldName = line[0]
            newName = line[1]
            origName = line[2]
            matches[newName] = oldName
            if origName in updated_names:
                updated_names[origName].append(oldName)
            else:
                updated_names[origName] = [oldName]

t.close()

# Search the new annotation for matching sequences, and rename the record using the old annotation names.

matching_seqs = {}
with open(new_annotations) as n:
    for record in SeqIO.parse(n, "fasta"):
        id = record.id
        if id in matches:
            newName = matches[id]
            record.id = newName
            record.name = newName
            record.description = newName
            record.seq = record.seq.lower()
            matching_seqs[newName] = record

n.close()

# Finally, replace the old sequences with the new sequences.

inhandle = open(old_annotations)
outhandle = open(final_annotations, "w")
count_merge = 0
count_split = 0

#print(updated_names.values())

for record in SeqIO.parse(inhandle, "fasta"):
    id = record.id
    if id in matching_seqs: 
        count_merge += SeqIO.write(matching_seqs[id], outhandle, "fasta")
    elif id in updated_names:
        for final_name in updated_names[id]:
            count_split += SeqIO.write(matching_seqs[final_name], outhandle, "fasta")
    else:
        SeqIO.write(record, outhandle, "fasta")

inhandle.close()
outhandle.close()
#print("%d sequences rewritten" % count)
