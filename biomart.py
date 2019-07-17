#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Use blast results to fill missing entries in the homology table.

@author: kaplanlab
"""

from pybiomart import Server
import pandas as pd
from os import path
import numpy as n
import subprocess 

tmp_folder = '/opt/current_project/tmp/eksikBlast/'

filenames = [path.join(tmp_folder, 'mapping_macaque.csv'), path.join(tmp_folder, 'mapping_human.csv'), 
             path.join(tmp_folder, 'mapping_chimp.csv')]

search_table = pd.concat([pd.read_csv(filename, header=None) for filename in filenames], axis=0)
search_table.columns = ['gene_id', 'protein_numbers']

def np_to_gene(np_id): 
    global search_table
    np_id = np_id.split('.')[0]
    
    row = search_table[search_table['protein_numbers'].str.contains(np_id)]['gene_id']
    
    if row.shape[0] == 0:
        return -1
    else:
        return row.iloc[0]


chimp_df = pd.read_csv(path.join(tmp_folder, 'results', 'Chimp.list'), sep='\t', header=None).iloc[:, :2]
chimp_df = chimp_df.drop_duplicates()

chimp_df.loc[:, 'human_gene_id'] = chimp_df.loc[:,0].apply(np_to_gene).astype(str)
chimp_df.loc[:, 'chimp_gene_id'] = chimp_df.loc[:, 1].apply(np_to_gene).astype(str)
chimp_df.loc[chimp_df[1] == 'XP_016789240.1', 'chimp_gene_id'] = '455432'
chimp_df.loc[chimp_df[1] == 'XP_016783511.1', 'chimp_gene_id'] = '453361'
chimp_df.loc[chimp_df[1] == 'XP_516182.3', 'chimp_gene_id'] = '460050'

chimp_processed_df = chimp_df.groupby(['human_gene_id'])['chimp_gene_id'].agg(lambda x: ','.join(set(x.dropna())))

update_rows = []

for human_gene_id, chimp_gene_id in chimp_processed_df.iteritems():
    update_rows.append("UPDATE homology SET chimp_gene_id='"+chimp_gene_id+"' WHERE "+
                       "human_gene_id='"+human_gene_id+"'" )   
    
chimp_update_sql = ';\n'.join(update_rows)

with open('/tmp/sanane_chimp.sql', 'w') as f:
    f.write(chimp_update_sql)


macaque_df = pd.read_csv(path.join(tmp_folder, 'results', 'Macaque.list'), sep='\t', header=None).iloc[:, :2]
macaque_df = macaque_df.drop_duplicates()

macaque_df.loc[:, 'human_gene_id'] = macaque_df.loc[:,0].apply(np_to_gene).astype(str)
macaque_df.loc[:, 'macaque_gene_id'] = macaque_df.loc[:, 1].apply(np_to_gene).astype(str)
macaque_df.loc[macaque_df[1] == 'XP_014965606.1', 'macaque_gene_id'] = '100430324'
macaque_df.loc[macaque_df[1] == 'XP_014996549.1', 'macaque_gene_id'] = '453361'
macaque_df.loc[macaque_df[1] == 'XP_014996287.1', 'macaque_gene_id'] = '692083'
macaque_df.loc[macaque_df[1] == 'XP_014996553.1', 'macaque_gene_id'] = '702071'

macaque_processed_df = macaque_df.groupby(['human_gene_id'])['macaque_gene_id'].agg(lambda x: ','.join(set(x.dropna())))

macaque_processed_df

update_rows = []

for human_gene_id, macaque_gene_id in macaque_processed_df.iteritems():
    update_rows.append("UPDATE homology SET macaque_gene_id='"+macaque_gene_id+"' WHERE "+
                       "human_gene_id='"+human_gene_id+"'" )   
    
macaque_update_sql = ';\n'.join(update_rows)

with open('/tmp/sanane.sql', 'w') as f:
    f.write(macaque_update_sql)
