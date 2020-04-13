
import glob2
import pandas as pd
from Bio import SeqIO
import pymysql
from tqdm import tqdm
import os
import re
import hashlib
import yaml
from pybiomart import Server
from pipeline.mysql_helper import MySqlHelper
from pipeline import connection_loader, config_loader
import logging
def load_mapping_collection():

    server = Server(host='http://www.ensembl.org')
    mart = server['ENSEMBL_MART_ENSEMBL']
    mouse_dataset = mart['mmusculus_gene_ensembl']
    human_dataset = mart['hsapiens_gene_ensembl']

    human_table = human_dataset.query(attributes=['ensembl_transcript_id', 'ensembl_peptide_id'])
    mouse_table = mouse_dataset.query(attributes=['ensembl_transcript_id', 'ensembl_peptide_id'])
    tables = pd.concat([human_table, mouse_table])
    tables.rename(columns={"Transcript stable ID": "transcript_id", "Protein stable ID": "protein_id"}, inplace=True)
    tables.dropna(how="any", inplace=True)

    mapping_collection = {}

    for index, row in tables.iterrows():
        prot_id = row["protein_id"]
        trans_id = row["transcript_id"]
        mapping_collection[prot_id] = trans_id

    return mapping_collection

def process_fasta_files(protein_path):
    species_list = ["Mus_musculus"]
    curation = {}
    species_dfs = []
    species_mapping_df = None

    mapping_collection = load_mapping_collection()

    # To collect the all fasta, faa, fa files in a dict
    for species in species_list:
        # species_mapping = []
        curation[species] = []
        
        for file in glob2.glob(os.path.join(protein_path, "*" + species + "*")):
            
            print(file)  
            for seq in tqdm(SeqIO.parse(file ,"fasta")):
                transcript_id = seq.id
                fasta = str(seq.seq).replace('\n','').replace("'", '')

                collector = {}
                collector["species"] = species
                if "|" in transcript_id:
                    transcript_id = transcript_id.split("|")[1]

                if transcript_id.startswith("ENSP") or transcript_id.startswith("ENSMUS"):
                    try:
                        transcript_id = mapping_collection[transcript_id.split(".")[0]]
                    except KeyError as e:
                        logging.error(f"ENSP-> ENST/ENSMUST mapping not found for {transcript_id}")
                        logging.error(f"Number of elements in collection {len(mapping_collection)}")
                        continue

                collector["db_id"] = transcript_id

                collector["seq"] = fasta
                # species_mapping.append(collector)
                curation[species].append(collector)

        species_mapping_df = pd.concat([species_mapping_df, pd.DataFrame(curation[species])])
        species_mapping_df.drop_duplicates(inplace=True)
        
    return species_mapping_df

def save_identifier_mapping_to_db(species_mapping_df, mysql_helper):
    for index, row in species_mapping_df.iterrows():
        fasta = row["seq"]
        db_id = row["db_id"]
        species = row["species"]
        mysql_helper.write_fasta_to_db(fasta, db_id, species.replace("_", " "))

if __name__ == '__main__':
    conn = connection_loader()
    mysql_helper = MySqlHelper(conn)
    config =  config_loader()

    species_mapping_df = process_fasta_files(os.path.join(config['staging_path'], 'proteins'))
    save_identifier_mapping_to_db(species_mapping_df, mysql_helper)