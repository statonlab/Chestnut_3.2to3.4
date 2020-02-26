#!/data/apps/python/3.2.1/bin/python
import re, sys, getopt
from Bio import SeqIO

# This script cleans up the remaining sequences by removing genes that got merged into other genes, isoforms that were lost, and any genes that were not called in the new pipeline.

input_annotations = "/staton/projects/chestnut/psudochro/analysis_081718_annotation/12k_RNA_Qrobur_081718/8_fixInternalStops/3_renameGenes/Castanea_mollissima_scaffolds_v3.4_HQcds_iter2.faa"

filter_out = "/staton/projects/chestnut/psudochro/analysis_081718_annotation/12k_RNA_Qrobur_081718/8_fixInternalStops/3_renameGenes/genes_to_filter.txt"

final_annotation = "/staton/projects/chestnut/psudochro/analysis_081718_annotation/12k_RNA_Qrobur_081718/8_fixInternalStops/3_renameGenes/Castanea_mollissima_scaffolds_v3.4_HQcds.faa"

# Create a list with every gene you want to filter out.

merged_list = []
with open(filter_out) as m:
    for line in m:
        merged_gene = line.rstrip()
        merged_list.append(merged_gene)

m.close()

# Use this list to filter out any sequences with a matching record ID:

inhandle = open(input_annotations)
outhandle = open(final_annotation, "w")
count = 0

for record in SeqIO.parse(inhandle, "fasta"):
    id = record.id
    if id not in merged_list:
        SeqIO.write(record, outhandle, "fasta")
    else:
        count += 1

inhandle.close
outhandle.close
print("%d genes were filtered out" % count)
