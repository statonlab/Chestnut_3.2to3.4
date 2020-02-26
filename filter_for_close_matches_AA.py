#!/data/apps/python/3.2.1/bin/python
import re, sys, getopt
from Bio import SeqIO

# This script will rename sequences from the updated annotation to match the original, and then replace the old predicted genes with the new ones.

input_tsv = "/staton/projects/chestnut/psudochro/analysis_081718_annotation/12k_RNA_Qrobur_081718/8_fixInternalStops/3_renameGenes/close_matches.txt"

old_annotations = "/staton/projects/chestnut/psudochro/analysis_081718_annotation/12k_RNA_Qrobur_081718/8_fixInternalStops/3_renameGenes/Castanea_mollissima_scaffolds_v3.2_HQcds.faa"

new_annotations = "/staton/projects/chestnut/psudochro/analysis_081718_annotation/12k_RNA_Qrobur_081718/8_fixInternalStops/3_renameGenes/augustus.hints.aa"

final_annotations = "/staton/projects/chestnut/psudochro/analysis_081718_annotation/12k_RNA_Qrobur_081718/8_fixInternalStops/3_renameGenes/Castanea_mollissima_scaffolds_v3.2_HQcds_iter1.faa"

matches = {}
with open(input_tsv) as t:
    for line in t:
        line = line.rstrip().split()
        oldName = line[0]
        newName = line[1]
        matches[newName] = oldName

t.close()

# Search new annotation for matching sequences, and rename the record using the old annotation names

matching_seqs = {}
with open(new_annotations) as n:
    for record in SeqIO.parse(n, "fasta"):
        id = record.id
        if id in matches:
            newName = matches[id]
            record.id = newName
            record.name = newName
            record.description = newName
            matching_seqs[newName] = record

n.close()

#print(matching_seqs)

# Finally, replace the old sequences with the new sequences.

inhandle = open(old_annotations)
outhandle = open(final_annotations, "w")
count = 0

for record in SeqIO.parse(inhandle, "fasta"):
    id = record.id
    if id in matching_seqs:
            count += SeqIO.write(matching_seqs[id], outhandle, "fasta")
    else:
        SeqIO.write(record, outhandle, "fasta")

inhandle.close()
outhandle.close()
print("%d sequences rewritten" % count)
