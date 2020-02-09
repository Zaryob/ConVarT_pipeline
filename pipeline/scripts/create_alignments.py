import sys
import os
from os import path
import pymysql
import re
import pandas as pd
import subprocess
import hashlib
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=pd.errors.ParserWarning)
warnings.simplefilter(action="ignore", category=pymysql.Warning)
import shutil
from sqlalchemy import create_engine
from multiprocessing import Pool
from tqdm import *
import time
import random
import string
import yaml
import argparse
import numpy as np 

class Options():
    
    def __init__(self):
        self.parse_arguments()      
        self.parse_config()
        self.validate_arguments()

    def parse_config(self):
        dir_path = os.path.dirname(os.path.realpath(__file__))

        with open(dir_path+'/../config.yml') as f:
            config = yaml.load(f)

        self.config = config
        
    def parse_arguments(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("--mode", choices=['construct_database','single_fasta'], required=True,
                            help="Mode of alignment: either `construct_database` or `single_fasta`")
        parser.add_argument("--truncate", type=bool, default=False,
                            help="Whether to truncate the DB or not if mode is `construct_database`")
        parser.add_argument("--input_path", type=str, nargs='+', default=[],
                            help='The input folder which contains fasta files which is used when mode is construct_database')      
                      
        parser.add_argument("--output_path", type=str, required=True, 
                        help="The output filename or path depending on the mode")

        args = parser.parse_args()
        self.args = args
        self.mode = args.mode
        self.truncate = args.truncate
        self.input_path = args.input_path
        self.output_path = args.output_path
    
    def validate_arguments(self):
        if self.mode == 'construct_database':
            for cur_path in self.input_path:
                if not os.path.exists(cur_path):
                    raise Exception('{} does not exist. `input_path` needs to contain fasta files that will be saved to DB'.format(cur_path))

def generate_tmp_file(size=10, chars=string.ascii_uppercase + string.digits): 
    filename = path.join('/tmp', ''.join(random.choice(chars) for x in range(size)))

    return filename

def catched_subprocess(command):
    try:
        return subprocess.check_output(command, shell=True)
    except Exception as e:
        print("Error", e, command)
        return ''.encode()

def get_pseudo_matches(filename, df_sequences):
    alignment_identity = df_sequences[['species']].copy()
    
    alignment_identity['scores'] = np.random.rand(alignment_identity.shape[0], 1)

    best_matches = alignment_identity[['scores', 'species']].groupby("species").idxmax() 

    return best_matches['scores']


def get_best_matches(filename, df_sequences):
    """
    Based on the fasta file containing genes from different species, the method outputs best aligned matches in other species for human gene  
    """

    dist_file = generate_tmp_file()
    catched_subprocess("cp " + filename + " " + dist_file+"fasta")

    alignment_output = catched_subprocess('clustalw -infile="'+dist_file+'fasta" -type=protein -outfile="'+
                        dist_file+'" -output=fasta -outorder=aligned').decode('utf-8')

    scores = re.findall("\(1:([0-9]+)\) Aligned. Score:  ([0-9]+)", alignment_output)
    scores.insert(0, ('1', 100))

    alignment_identity = df_sequences[['species']].copy()
    
    alignment_identity['scores'] = [int(score_tuple[1]) for score_tuple in scores]

    best_matches = alignment_identity[['scores', 'species']].groupby("species").idxmax() 

    return best_matches['scores']

def process_fasta_file(filename):
    """
    Returns a data frame containing sequences in the fasta file and two columns: species and sequence, and the indices of rows denote
    np id.
    """
    sequences_list = []

    with open(filename, "r") as f:
        all_seqs = ('\n'+''.join(f.readlines()).replace(r"(>){1,}",">")).split("\n>")
        
        for i, cur_seq in enumerate(all_seqs):
            if cur_seq == "":
                continue

            cur_seq = '\n>'+cur_seq
            species_name = re.search(r"\[([^\[]+?)\]\n", cur_seq).group(1)
            protein_id = re.search(r">([^ ]+)", cur_seq).group(1)

            sequences_list.append({'species': species_name, 'np_id': protein_id, 'sequence': cur_seq})

    df_sequences = pd.DataFrame(sequences_list)
    if(df_sequences.shape[0] != 0):
        df_sequences = df_sequences.set_index('np_id')
    
    return df_sequences

def get_combinations(filename, df_sequences):
    
    #best_matches = get_pseudo_matches(filename, df_sequences)
    
    
    all_possible_combinations = []
    species_list = df_sequences['species'].unique().tolist()
    root_match = pd.Series(data=[str(df_sequences[df_sequences['species']==species_list[0]].index[0])], index=[species_list[0]])
    stack = [root_match]
    while len(stack) > 0:
        current_combination = stack.pop().copy()
        level = current_combination.shape[0]
        if level == len(species_list):
            all_possible_combinations.append(current_combination)
            continue
        
        for np_id in df_sequences[df_sequences['species'] == species_list[level]].index.tolist():
            current_combination.loc[species_list[level]] = str(np_id)
            stack.append(current_combination.copy())

    return all_possible_combinations

def create_alignment(input_file, df_sequences=None):
    output_file = generate_tmp_file()
    catched_subprocess('clustalw -infile='+input_file+' -type=protein -outfile='+output_file+' -output=fasta -outorder=aligned')

    with open(output_file, 'r') as file:
        fasta = ''.join(file.readlines())

    if type(df_sequences).__name__ != 'NoneType':
        for np_id in df_sequences.index.tolist():
            fasta = fasta.replace(np_id, np_id + ' ['+df_sequences.loc[np_id, 'species']+']')
        
        fasta = fasta.replace(r"(>){1,}", ">")
        fasta_list = ('\n'+fasta).split('\n>')
        fasta_list = list(filter(lambda x: x != '', fasta_list))

        human_gene_ind = fasta_list.index(list(filter(lambda x: 'Homo sapiens' in x, fasta_list))[0])

        human_sequence = fasta_list[human_gene_ind]
        del fasta_list[human_gene_ind]
        fasta_list.insert(0, human_sequence)

        fasta_txt = '>' + '\n>'.join(fasta_list)
    else:
        fasta_txt = fasta

    return fasta_txt

def align_combination(output_dir, combination, df_sequences):

    ordered_species = ["Homo sapiens", "Pan troglodytes", "Macaca mulatta", "Mus musculus", "Rattus norvegicus", "Danio rerio", 
                             "Xenopus tropicalis", "Caenorhabditis elegans", "Drosophila melanogaster"]

    fasta_content = ""

    for species_name in ordered_species:
        if df_sequences[df_sequences['species']==species_name].shape[0] == 0:
            continue

        fasta_content += df_sequences.loc[combination.loc[species_name], 'sequence']

        if species_name != ordered_species[-1]:
            fasta_content += "\n\n"
    
    new_filename = generate_tmp_file()
    
    processed_fasta = path.join(output_dir, 'processed_fasta', new_filename)

    with open(processed_fasta, "w") as f:
        f.write(fasta_content)

    return create_alignment(processed_fasta, df_sequences.loc[combination.tolist(), :])



def write_msa_to_db(fasta, combination):
    global cur, con

    cur.execute("INSERT INTO msa (fasta, alignment_method) VALUES(%s, %s)", (fasta, 'clustalw'))
    
    msa_id = cur.lastrowid

    for np_id in combination.values:
        np_id_without_version = np_id.split('.')[0]

        cur.execute("SELECT convart_gene_id FROM convart_gene_to_db WHERE db_id = %s LIMIT 1", 
                                        (np_id_without_version))
        row = cur.fetchone()
        convart_gene_id = row[0]
        con.commit()

        cur.execute("INSERT INTO msa_gene (msa_id, convart_gene_id) VALUES(%s, %s)", (msa_id, convart_gene_id))
    
    con.commit()

    return msa_id

def write_fasta_to_db(fasta, db_id, species):
    global con, cur
    
    db_id_version = None
    if '.' in db_id:
        db_id, db_id_version = db_id.split('.')[0], db_id.split('.')[1]

    if 'NP' in db_id:
        db = 'NP'
    elif 'ENST' in db_id:
        db = 'ENST'
    else:
        db = 'OTHER'

    fasta = re.sub(r'>(.*?)\n', '', fasta).replace('\n','').replace("'", '')
    seq_hash = hashlib.md5((fasta+species).encode('utf-8')).hexdigest()

    cur.execute("SELECT id FROM convart_gene WHERE hash=%s and species_id=%s", (seq_hash, species))
    
    if cur.rowcount == 0:
        cur.execute("INSERT IGNORE INTO convart_gene (sequence, species_id, hash) VALUES(%s, %s, %s)", 
                            (fasta, species, seq_hash))
        
        convart_gene_id = con.insert_id()
        con.commit()
    else:
        row = cur.fetchone()
        convart_gene_id = row[0]
    cur.execute("INSERT IGNORE INTO convart_gene_to_db (convart_gene_id, db, db_id, db_id_version) VALUES(%s, %s, %s, %s)", 
                    (convart_gene_id, db, db_id, db_id_version))

    return convart_gene_id

def write_best_combination_to_db(msa_id, convart_gene_id):
    global con, cur
    
    cur.execute("INSERT IGNORE INTO msa_best_combination (msa_id, convart_gene_id) VALUES(%s, %s)", 
                            (msa_id, convart_gene_id))
    con.commit()


def save_to_db(df_sequences, combinations, msa_results):
    for gene_id, row in df_sequences.iterrows():
        
        convart_gene_id = write_fasta_to_db(row['sequence'], gene_id, row['species'])
        if row['species'] == 'Homo sapiens':
            human_convart_gene_id = convart_gene_id
        
    for ind, combination in enumerate(combinations):
        fasta = msa_results[ind]
        msa_id = write_msa_to_db(fasta, combination)
        #if gene_id == df_sequences.index[0]:
        #    write_best_combination_to_db(msa_id, human_convart_gene_id)

## Creating fasta files that contains at most one gene for each species
def align_fasta_file(file):
    try:
        #source = re.search(r"_(.*?)\.fasta", file).group(1)
        
        filename = path.join(input_dir, file)
        df_sequences = process_fasta_file(filename)
        if df_sequences.shape[0] <= 1:
            return False

        if len(df_sequences['species'].unique()) == df_sequences.shape[0]:
            combination = pd.Series(df_sequences.index, index=df_sequences['species'])
            combinations = list({combination.loc['Homo sapiens']: combination}.values())
        else:
            combinations = get_combinations(filename, df_sequences)        
        
        msa_results = {}

        for ind, combination in enumerate(combinations):
            fasta = align_combination(output_dir, combination, df_sequences)
            msa_results[ind] = fasta
        
        return df_sequences, combinations, msa_results
                        
    except Exception as e:
        print(e, file)

def align_and_save_parallel(func, args, n_processes = 1):
    p = Pool(n_processes)

    with tqdm(total = len(args)) as pbar:
        for res in tqdm(p.imap_unordered(func, args, chunksize=n_processes)):
            pbar.update()
            if type(res) == bool or res == None:
                continue

            try:
                df_sequences, combinations, msa_results = res
                save_to_db(df_sequences, combinations, msa_results)
            except Exception as e:
                print("Error during insertion to db", sys.exc_info()[0], e)

    pbar.close()
    p.close()
    p.join()

if __name__ == '__main__':
    pipeline_options = Options()
    
    con = pymysql.connect(host=pipeline_options.config['MYSQL_HOST'], user=pipeline_options.config['MYSQL_USER'], passwd=pipeline_options.config['MYSQL_PASSWD'], 
                      db=pipeline_options.config['MYSQL_DB'])

    cur = con.cursor()
        
    mode = pipeline_options.mode

    if mode == 'construct_database':
        input_dirs = pipeline_options.input_path
        output_dir = pipeline_options.output_path

        if pipeline_options.truncate:
            
            print("Truncating MSA tables")

            cur.execute('TRUNCATE TABLE msa')
            cur.execute('TRUNCATE TABLE msa_gene')
            cur.execute('TRUNCATE TABLE msa_best_combination')
            cur.execute('TRUNCATE TABLE convart_gene')
            cur.execute('TRUNCATE TABLE convart_gene_to_db')
            con.commit()
        else:
            print("Inserting the MSA tables without truncating them")

        #Directory which stores the processed fasta files which contains at most one gene for each species 
        processed_fasta_output_dir = path.join(output_dir, 'processed_fasta')
        
        if  os.path.exists(processed_fasta_output_dir):
            shutil.rmtree(processed_fasta_output_dir)
        
        os.makedirs(processed_fasta_output_dir)
        
        for input_dir in input_dirs:
            raw_fasta_files = os.listdir(input_dir)
            results = align_and_save_parallel(align_fasta_file, raw_fasta_files, 27)

    elif mode == 'single_fasta':
        filename = pipeline.input_path[0]
        if not os.path.exists(filename):
            print("File does not exist")


        print(create_alignment(filename))
