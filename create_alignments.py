import sys
import os
from os import path
import pymysql
import re
import pandas as pd
import subprocess
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

def generate_tmp_file(size=10, chars=string.ascii_uppercase + string.digits): 
    filename = path.join('/tmp', ''.join(random.choice(chars) for x in range(size)))

    return filename

def catched_subprocess(command):
    try:
        return subprocess.check_output(command, shell=True)
    except Exception as e:
        print("Error", e, command)
        return ''.encode()

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
    
    best_matches = get_best_matches(filename, df_sequences)
    
    combinations = {}
    combinations[df_sequences.index[0]] = best_matches
    #print("Warning! we only consider best matching combination")
    
    return combinations

    for species_name in df_sequences['species'].unique().tolist():
        for np_id in df_sequences[df_sequences['species'] == species_name].index.tolist():
            if np_id == best_matches.loc[species_name]:
                continue

            fasta_combination = best_matches.copy()
            fasta_combination.loc[species_name] = np_id
            combinations[np_id] = fasta_combination

    return combinations

def create_alignment(input_file, df_sequences=None):
    output_file = generate_tmp_file()
    catched_subprocess('clustalw -infile='+input_file+' -type=protein -outfile='+output_file+' -output=fasta -outorder=aligned')

    with open(output_file, 'r') as file:
        fasta = ''.join(file.readlines())

    if type(df_sequences) != None:
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

con = pymysql.connect(host='127.0.0.1', unix_socket='/opt/lampp/var/mysql/mysql.sock', 
            user='root', passwd='', db='current_project')
cur = con.cursor()

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

    fasta = re.sub(r'>(.*?)\n', '', fasta).replace('\n','')
    cur.execute("SELECT id FROM convart_gene WHERE sequence = %s AND species_id= %s", (fasta, species))
    
    if cur.rowcount == 0:
        cur.execute("INSERT IGNORE INTO convart_gene (sequence, species_id) VALUES(%s, %s)", 
                            (fasta, species))
        
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
    
    cur.execute("INSERT INTO msa_best_combination (msa_id, convart_gene_id) VALUES(%s, %s)", 
                            (msa_id, convart_gene_id))
    con.commit()


def save_to_db(df_sequences, combinations, msa_results):
    for gene_id, row in df_sequences.iterrows():
        convart_gene_id = write_fasta_to_db(row['sequence'], gene_id, row['species'])
        if row['species'] == 'Homo sapiens':
            human_convart_gene_id = convart_gene_id

    for gene_id, combination in combinations.items():
        fasta = msa_results[gene_id]
        msa_id = write_msa_to_db(fasta, combination)
        if gene_id == df_sequences.index[0]:
            write_best_combination_to_db(msa_id, human_convart_gene_id)

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
            combinations = {combination.loc['Homo sapiens']: combination}
        else:
            combinations = get_combinations(filename, df_sequences)        
        
        msa_results = {}

        for protein_id, combination in combinations.items():
            fasta = align_combination(output_dir, combination, df_sequences)
            msa_results[protein_id] = fasta
        
        return df_sequences, combinations, msa_results
                        
    except Exception as e:
        print(e, file)

def align_and_save_parallel(func, args, n_processes = 30):
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
    if len(sys.argv) < 2:
        print("The correct format is python3 create_alignment.py <command (i.e. construct_database single_fasta)> <filename if command is single_fasta>")

    mode = sys.argv[1]

    if mode == 'construct_database':
        input_dirs = ["/opt/current_project/results/seqs_with_homology"]
        output_dir = "/opt/current_project/results/alignments_new3"

        #Directory which stores the processed fasta files which contains at most one gene for each species 
        processed_fasta_output_dir = path.join(output_dir, 'processed_fasta')
        
        if  os.path.exists(processed_fasta_output_dir):
            shutil.rmtree(processed_fasta_output_dir)
        
        os.makedirs(processed_fasta_output_dir)

        cur.execute('TRUNCATE TABLE msa')
        cur.execute('TRUNCATE TABLE msa_gene')
        cur.execute('TRUNCATE TABLE msa_best_combination')
        cur.execute('TRUNCATE TABLE convart_gene')
        cur.execute('TRUNCATE TABLE convart_gene_to_db')
        con.commit()
        
        for input_dir in input_dirs:
            raw_fasta_files = os.listdir(input_dir)
            results = align_and_save_parallel(align_fasta_file, raw_fasta_files, 30)

    elif mode == 'single_fasta':
        if len(sys.argv) != 3 or not os.path.exists(sys.argv[2]):
            print("File does not exist")

        filename = sys.argv[2]

        print(create_alignment(filename))
