#!/data/apps/python/3.2.1/bin/python
import re, sys, getopt
from Bio import SeqIO

# This script will filter out genes that merged into other genes, isoforms that were lost, and any gene that didn't get called in the new pipeline.

input_annotations = "/staton/projects/chestnut/psudochro/analysis_081718_annotation/12k_RNA_Qrobur_081718/8_fixInternalStops/3_renameGenes/Castanea_mollissima_scaffolds_v3.4_HQ_iter2.gff"

merged_out = "/staton/projects/chestnut/psudochro/analysis_081718_annotation/12k_RNA_Qrobur_081718/8_fixInternalStops/3_renameGenes/merged_genes.txt"

lost_isoforms = "/staton/projects/chestnut/psudochro/analysis_081718_annotation/12k_RNA_Qrobur_081718/8_fixInternalStops/3_renameGenes/lost_isoforms.txt"

junk_genes = "/staton/projects/chestnut/psudochro/analysis_081718_annotation/12k_RNA_Qrobur_081718/8_fixInternalStops/3_renameGenes/junk_genes.txt"

final_annotation = "/staton/projects/chestnut/psudochro/analysis_081718_annotation/12k_RNA_Qrobur_081718/8_fixInternalStops/3_renameGenes/Castanea_mollissima_scaffolds_v3.4_HQ.gff"

# Since you need to approach each of the lists differently, I've separated them into separate lists.

# For example, you will want the parent names of the merged out genes and the junk genes, but not the lost isoforms.

merged_list = []
with open(merged_out) as m:
    for line in m:
        merged_gene = line.rstrip()
        merged_parent = re.sub(".t\d","",merged_gene)
        merged_list.append(merged_parent)

m.close()

lost_list = []
with open(lost_isoforms) as i:
    for line in i:
        isoform = line.rstrip()
        lost_list.append(isoform)

junk_list = []
with open(junk_genes) as j:
    for line in j:
        junk_gene = line.rstrip()
        junk_parent = re.sub(".t\d","",junk_gene)
        if junk_parent not in junk_list:
            junk_list.append(junk_parent)

# With your three lists, we're going to filter all of these genes out of the GFF file.

inhandle = open(input_annotations)
outhandle = open(final_annotation, "w")

for anno in inhandle:
    line = anno.rstrip().split()
    if line[2] == "gene":
        gene = line[8].split(";")
        gene = re.sub(r"ID=","",gene[0])
        if (gene not in junk_list) and (gene not in merged_list):
            outhandle.write(anno)
    else:
        gene = line[8].split(";")
        feature = re.sub(r"ID=","",gene[0])
        parent = re.sub(r"Parent=","",gene[1])
        if line[2] == "mRNA":
            if (feature not in lost_list) and (parent not in junk_list) and (parent not in merged_list):
                outhandle.write(anno)
        elif line[2] == "CDS":
            feature = re.sub(r".CDS\d+","",feature)
            parent = re.sub(r".t\d+","",parent)
            if (feature not in lost_list) and (parent not in junk_list) and (parent not in merged_list):
                outhandle.write(anno)
        else:
            feature = re.sub(r".exon\d+","",feature)
            parent = re.sub(r".t\d+","",parent)
            if (feature not in lost_list) and (parent not in junk_list) and (parent not in merged_list):
                outhandle.write(anno)

inhandle.close()
outhandle.close()
