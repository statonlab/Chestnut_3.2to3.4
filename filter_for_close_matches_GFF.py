#!/data/apps/python/3.2.1/bin/python
import re, sys, getopt
from Bio import SeqIO

# This script will rename sequences from the updated annotation to match the original, and then replace the old predicted genes with the new ones.

input_tsv = "/staton/projects/chestnut/psudochro/analysis_081718_annotation/12k_RNA_Qrobur_081718/8_fixInternalStops/3_renameGenes/close_matches.txt"

old_annotations = "/staton/projects/chestnut/psudochro/analysis_081718_annotation/12k_RNA_Qrobur_081718/8_fixInternalStops/3_renameGenes/Castanea_mollissima_scaffolds_v3.2_HQ.gff"

new_annotations = "/staton/projects/chestnut/psudochro/analysis_081718_annotation/12k_RNA_Qrobur_081718/8_fixInternalStops/3_renameGenes/augustus.hints.gff3"

final_annotations = "/staton/projects/chestnut/psudochro/analysis_081718_annotation/12k_RNA_Qrobur_081718/8_fixInternalStops/3_renameGenes/Castanea_mollissima_scaffolds_v3.4_HQ_iter.gff"

matches = {}
parentMatches = {}
with open(input_tsv) as t:
    for line in t:
        line = line.rstrip().split()
        oldName = line[0]
        oldParent = re.sub(".t\d","",oldName)
        newName = line[1]
        newParent = re.sub(".t\d","",newName)
        matches[newName] = oldName 
        parentMatches[newParent] = oldParent

t.close()

gffFix = {}

with open(new_annotations) as n:
    for line in n:
        line = line.rstrip().split()
        if line[2] == "gene":
            gene_name = line[8].split(";")
            gene_name = re.sub(r"ID=","",gene_name[0])
            if gene_name in parentMatches:
                origName = parentMatches[gene_name]
                line[8] = line[8].replace(gene_name, origName)
                gffFix[origName] = ["\t".join(line)]
        elif line[2] == "mRNA":
            gene_name = line[8].split(";")
            gene_name = re.sub(r"ID=","",gene_name[0])
            if gene_name in matches:
                origName = matches[gene_name]
                parentGene = re.sub(".t\d+","",gene_name)
                parentName = re.sub(".t\d+","",origName)
                line[8] = line[8].replace(gene_name, origName)
                line[8] = line[8].replace(parentGene, parentName)
                gffFix[parentName].append("\t".join(line))
        else:
            gene_name = line[8].split(";")
            gene_name = re.sub(r"Parent=","",gene_name[1])
            if gene_name in matches:
                origName = matches[gene_name]
                parentGene = re.sub(".t\d+","",gene_name)
                parentName = re.sub(".t\d+","",origName)
                line[8] = line[8].replace(gene_name, origName)
                line[8] = line[8].replace(parentGene, parentName)
                if line[2] in ["CDS", "exon"]:
                    gffFix[parentName].append("\t".join(line))

n.close()

inhandle = open(old_annotations)
outhandle = open(final_annotations, "w")

for anno in inhandle:
    if "#" not in anno:
        line = anno.rstrip().split()
        if line[2] == "gene":
            gene = line[8].split(";")
            gene = re.sub(r"ID=","",gene[0])
            if gene in gffFix:
                for l in gffFix[gene]: 
                    outhandle.write(l + "\n")
            else:
                outhandle.write(anno)
        else:
            gene = line[8].split(";")
            gene = re.sub(r"Parent=","",gene[1])
            if gene not in gffFix:
                if line[2] == "mRNA":
                    outhandle.write(anno)
                else:
                   gene = re.sub(".t\d+","",gene)
                   if gene not in gffFix:
                       outhandle.write(anno)

inhandle.close()
outhandle.close()
