#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 13:07:36 2020

@author: annieforster
"""

import pandas as pd
import sys
import os

# script to identify the genomic region type (e.g. genic, intergenic, 3' UTR etc.) of control regions generated using control_seq_random_bed.py
# first run bedtools intersect on control regions generated 
# e.g bedtools intersect -a Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.12.gff3 -b <control_seg_bed_file> -bed > <control_seq_intersection_bed_file>

control_seq = sys.argv[1] # control regions generated by control_seq_random_bed.py
intersect_file = sys.argv[2] # intersection of control regions and Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.12.gff3 using bedtools intersect

df_original = pd.read_csv(intersect_file, sep='\t', header=None)
df = df_original[df_original[2] != "chromosome"]
df = df[df[2] != "CDS"]
df = df[df[2] != "mRNA"]
df = df.drop_duplicates(subset=[0, 1, 2, 3, 4], inplace=False)

sequences = open(control_seq, "r")

number = sys.argv[3]

counter = 0
chromosome = str(control_seq).split("_")[0]
filename = chromosome + "_control_sequence_categories_" + str(number) + ".txt"

output = open(filename, "w")

counter = 0
chromosome = str(control_seq).split("_")[0]
filename = chromosome + "_control_sequence_MANUAL_CHECK_" + str(number) + ".txt"

manual_check = open(filename, "w")

output.write("Start\tEnd\tGO\n")

def original():
    original_locus = df_original.loc[(df_original[4] >= start_pos) & (df_original[3] <= start_pos) | (df_original[4] >= end_pos) & (df_original[3] <= end_pos) | (df[4] <= end_pos) & (df[3] >= start_pos)]
    return original_locus

def distinct(x):
    global previous
    if previous != x:
        if not count == 0:
            output.write("/")
        output.write(str(x))
    previous = x

n = 1
for sequence in sequences.readlines():
    end = 0
    positions = sequence.rstrip("\n").split("\t")
    start_pos = int(positions[1]) + 1
    end_pos = int(positions[2]) 
    output.write(str(start_pos) + "\t" + str(end_pos) + "\t")
    locus = df.loc[(df[4] >= start_pos) & (df[3] <= start_pos) | (df[4] >= end_pos) & (df[3] <= end_pos) | (df[4] <= end_pos) & (df[3] >= start_pos)]
    check = locus.loc[(locus[4] == locus[3])]
    if not check.empty:
        output.write(str(n) + ". Manual Check\n")
        manual_check.write(str(n) + ".\n" + str(original()) + "\n")         
        n += 1
        continue
    if locus.empty:
        output.write("intergenic\n")
        continue
    original_locus = original()
    len_mRNA = len(original_locus.loc[original_locus[2] == "mRNA"])
    if (len_mRNA >= 2):
        exon_locus = original_locus.loc[original_locus[2] == "exon"]
        len_exon = len(exon_locus)
        if (exon_locus[4].nunique() >= 2) or (exon_locus[3].nunique() >= 2):

            output.write(str(n) + ". Manual Check - discordant overlapping mRNAs\n")
            manual_check.write(str(n) + ".\n" + str(original()) + "\n")         
            n += 1
            continue
    gene_locus = locus.loc[(locus[2] == "gene")]
    if gene_locus.empty:
        not_gene_locus = locus.loc[(locus[2] == "pre_miRNA") | (locus[2] == "tRNA")]
        try: 
            gene_type = not_gene_locus[4].item()
            if gene_type != end_pos:
                end = 1
        except ValueError:
            output.write(str(n) + ".Manual Check: not a gene and overlapping\n")
            manual_check.write(str(n) + ".\n" + str(original()) + "\n")         
            n += 1
            continue
    else:
        if len(gene_locus) >= 2:
            output.write(str(n) + ".Manual Check: more than one overlapping gene\n")
            manual_check.write(str(n) + ".\n" + str(original()) + "\n")         
            n += 1
            continue
        if gene_locus[4].iloc[0] != end_pos:
            end = 1
    sections1 = locus.drop_duplicates(subset=[3], inplace=False)
    sections2 = locus.drop_duplicates(subset=[4], inplace=False)
    start_points = sections1[3].tolist()
    end_points = sections2[4].tolist()
    points = start_points + end_points
    final_points = [int(i) for i in points]
    final_points.sort()
    if final_points[0] != start_pos:
        output.write("intergenic/")
    locus = locus[locus[2] != "gene"]
    locus = locus[locus[2] != "ncRNA_gene"]
    nums = len(final_points)
    count = 0 
    previous = ""      
    for elem,next_elem in zip(final_points, final_points[1:]+[final_points[0]]):
        if elem == final_points[-1]:
            break
#        new_locus = locus.loc[(locus[3] <= next_elem - 1) & (locus[3] >= elem)]
        new_locus = locus.loc[(locus[3] == elem) & (locus[4] != next_elem) | (locus[3] != elem) & (locus[4] == next_elem) | (locus[3] == elem) & (locus[4] == next_elem)]
        while True:
            if new_locus.empty:
                distinct("intron")
                break
            elif len(new_locus) == 1:
                distinct(str(new_locus[2].item()))
                break
            elif len(new_locus) >= 2:
                if new_locus[2].nunique() == 1:
                    GO = new_locus[2].iloc[0]
                    distinct(str(GO))
                    break
                filter_trial = new_locus[new_locus[2] != "exon"]
                filter_trial = filter_trial[filter_trial[2] != "three_prime_UTR"]
                filter_trial = filter_trial[filter_trial[2] != "five_prime_UTR"]
                if filter_trial.empty:
                    new_locus = new_locus[new_locus[2] != "exon"]
                else:
                        if len(filter_trial) == 1:
                            new_locus = filter_trial
                        else:
                            output.write(str(n) + "Manual Check. rare overlapping variants\n")
                            manual_check.write(str(n) + ".\n" + str(original()) + "\n")         
                            n += 1
                            break                       
            else:
                output.write(str(n) + ". Manual Check")
                manual_check.write(str(n) + ".\n" + str(original()) + "\n")         
                n += 1
                break
        count += 1
  
    if end == 1:
        output.write("/intergenic")
        end = 0
    output.write("\n")
        
output.close()
    