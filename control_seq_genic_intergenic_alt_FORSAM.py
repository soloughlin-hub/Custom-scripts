#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 13:48:30 2020

@author: annieforster
"""

# generates a set of control regions for the invariant regions in the Anopheles gambiae genome
# matching for whether the invariant region fell in the genic or intergenic region of the genome.

from random import random 
import pandas as pd
import sys
import os

ACEs_in = sys.argv[1] # tab delimited bed file of invariant region IDs for a single chromosome, start postion, end position and whether was in genic or intergenic region on chromosome
# e.g. ID1\t1000\t2000\tgenic
AGAPs = sys.argv[2] # tab delimited  file of AGAP gene IDs for a single chromosome, start and end positions. GFF3 file f or single chromosome arm.
intergenics = sys.argv[3] # tab delimited file of intergenic regions on chromosome
# e.g. 2L\t1000\t2000

chromosome = ACEs_in.split("_")[0]

counter = 0
filename = str(chromosome) + "_Control_Regions_2_Positions_alt_{}.bed"
while os.path.isfile(filename.format(counter)):
    counter += 1
if counter == 0:
    filename = str(chromosome) + "_Control_Regions_2_Positions_alt_0.bed"
else:
    filename = filename.format(counter)
output = open(filename, "w")

length_d = dict([
    ('2L', 49364325),
    ('2R', 61545105),
    ('3L', 41963435), 
    ('3R',  53200684),
    ('X', 24393108)])

length = length_d[chromosome]

df_int = pd.read_csv(intergenics, sep='\t', header=None)
df_int['Length'] = (df_int[2] - df_int[1]) + 1
df_gen = pd.read_csv(AGAPs, sep='\t', header=None)
df_gen['Length'] = (df_gen[2] - df_gen[1]) + 1

ACEs = open(ACEs_in, "r")

df2 = pd.DataFrame(columns=[0, 1])

for line in ACEs.readlines():
    parts = line.split("\t")
    gen_type = parts[3].rstrip("\n")
    start = parts[1]
    end = parts[2]
    ACE_length = (int(end) - int(start)) + 1
    if gen_type == "genic":
        df = df_gen
    else:
        df = df_int
    while True:
        total_row = len(df.index)
        ran = random()
        chunk = int(total_row * ran) - 1
        locus = df.iloc[chunk, :]
        loc_len = locus['Length'].item()
        if loc_len < ACE_length:
            pass
        else:
            ran2 = random()
            restricted = (loc_len - ACE_length) + 1
            select = int(ran2 * restricted)
            position = locus[1].item() + select
            control_end = position + (ACE_length - 1)
            locus_check = df2.loc[((df2[0] >= position) & (df2[0] <= control_end)) | ((df2[1] >= position) & (df2[1] <= control_end)) | ((df2[0] <= position) & (df2[1] >= control_end)) | ((df2[0] >= position) & (df2[1] <= control_end))]
            if locus_check.empty:
                if control_end <= locus[2].item():
                    output.write(chromosome + "\t" + str(int(position)) + "\t" + str(int(control_end)) + "\t" + gen_type + "\n")
                    df2.loc[len(df2)] = [int(position),int(control_end)]
                    break
                else:
                    print("error")
        
ACEs.close()
output.close()
  