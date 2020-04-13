#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import glob2

files = glob2.glob("db/tubulin_mut/" + "*.csv")
dfs = [pd.read_csv(file, quotechar='"', usecols=["Organism", "Tubulin Gene (Isotype)", "AA3", "Phenotype"]) for file in files]
dfs = [df.dropna(how="all") for df in dfs]

for df in dfs:
    df["position"] = df["AA3"].str.extract('([0-9][0-9]*)', expand=True)
    df["transcript"] = df["Tubulin Gene (Isotype)"].apply(lambda x: x.split(" ")[0])
    df = df.rename(columns={'AA3':'Aa_change', 'Tubulin Gene (Isotype)': "gene_tag"}, inplace=True)

for i in range(len(dfs)):
    dfs[i] = dfs[i][dfs[i]["Organism"].isin(["H. sapiens", "M. musculus", "C. elegans"])]

final = dfs[0].append(dfs[1].append(dfs[2]))
final.to_csv("tubulin_mutations_3species.csv", sep="\t", index=False)

