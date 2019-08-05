#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 20:29:33 2019

@author: kaplanlab
"""
import pandas as pd
from matplotlib import pyplot as plt
import pymysql

con = pymysql.connect(host='127.0.0.1', unix_socket='/opt/lampp/var/mysql/mysql.sock', user='root', passwd='', db='current_project')
cur = con.cursor()

df = pd.read_sql("SELECT clinvar_id, gene_id, name FROM clinvar WHERE "+
                 "gene_id!=0", con)

def extract_nm_id_from_variant_name(name):
    nm_id_with_version = name.split('(')[0]
    
    nm_id_without_version = nm_id_with_version.split('.')[0]
    
    return nm_id_without_version

df.loc[:, 'nm_id'] = df['name'].apply(extract_nm_id_from_variant_name)

sql = "UPDATE clinvar SET nm_id='{}' WHERE clinvar_id={};"
queries = df.apply(lambda x: sql.format(x['nm_id'], x['clinvar_id']), axis=1)

for query in queries:
    cur.execute(query)
