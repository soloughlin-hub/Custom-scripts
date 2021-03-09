#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 16:35:57 2019

@author: annieforster
"""

# generates a set of control sequences for the invariant regions in Anopheles gambiae genome matching 
# the set in number and invariant region length

from random import random 
import pandas as pd
import sys
import os

IR_file_in = sys.argv[1] # tab delimited text file with invariant region ID and invariant region length
# e.g. ID1\t100

chromosome = str(IR_file_in).split("_")[0]

counter = 0
filename = str(chromosome) + "_Control_Regions_Positions_{}.bed"
while os.path.isfile(filename.format(counter)):
    counter += 1
if counter == 0:
    filename = str(chromosome) + "_Control_Regions_Positions_0.bed"
else:
    filename = filename.format(counter)
output = open(filename, "w")

index_file_name = "/chr" + chromosome + "_control_seq_index.txt"
directory = os.getcwd()
index_file = directory + index_file_name
df = pd.read_csv(index_file, sep='\t', header=0)
df['Length'] = df['Length'].astype(int)
df['Start_Pos'] = df['Start_Pos'].astype(int)
df['End_Pos'] = df['End_Pos'].astype(int)
df['Start_Range'] = df['Start_Range'].astype(float)
df['End_Range'] = df['End_Range'].astype(float)

IR_file = open(IR_file_in, "r")

df2 = pd.DataFrame(columns=[0, 1])

for line in IR_file.readlines():
    while True:
        length = (line.rstrip("\n")).split("\t")[1]
        ran = random()
        locus = (df.loc[(df['End_Range'] >= ran) & (df['Start_Range'] <= ran)])
        into = round((locus['Length'].item())*ran, 0)
        start_pos = into + locus['Start_Pos'].item()
        length_selection = int(start_pos) + int(length) - 1
        if length_selection < locus['End_Pos'].item():
            locus_check = df2.loc[(df2[0] >= start_pos) & (df2[0] <= length_selection) | (df2[1] >= start_pos) & (df2[1] <= length_selection) | (df2[0] <= start_pos) & (df2[1] >= length_selection) | (df2[0] >= start_pos) & (df2[1] <= length_selection)]
            if locus_check.empty:
                output.write(chromosome + "\t" + str(int(start_pos)) + "\t" + str(length_selection) + "\n")
                df2.loc[len(df2)] = [int(start_pos),int(length_selection)]
                break
            else:
                pass
        else:
            pass 

output.close()
IR_file.close()
    