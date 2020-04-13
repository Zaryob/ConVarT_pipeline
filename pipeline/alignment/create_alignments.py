import multiprocessing
multiprocessing.set_start_method('spawn', True)

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

from tqdm import *
import time
import random
import string
import argparse
import numpy as np 
from pipeline.mysql_helper import MySqlHelper
from pipeline.parallel_helper import run_parallel
from pipeline import (
    config_loader,
    catched_subprocess,
    connection_loader
)

class Options():
    
    def __init__(self):
        self.parse_arguments()      
        self.parse_config()
        self.validate_arguments()

    def parse_config(self):
        self.config = config_loader()

    def parse_arguments(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("--mode", choices=['construct_database','single_fasta'], required=True,
                            help="Mode of alignment: either `construct_database` or `single_fasta`")
        parser.add_argument("--database", type=str, required=True,
                            help="Database name")
        parser.add_argument("--truncate", type=bool, default=False,
                            help="Whether to truncate the DB or not if mode is `construct_database`")
        parser.add_argument("--input_path", type=str, nargs='+', default=[],
                            help='The input folder which contains fasta files which is used when mode is construct_database')      
        parser.add_argument("--only_best_match", type=bool, default=False,
                            help="If set true, creates MSAs for only the best match. Default is False")
        parser.add_argument("--output_path", type=str, required=True, 
                        help="The output filename or path depending on the mode")

        args = parser.parse_args()
        self.args = args
        self.mode = args.mode
        self.truncate = args.truncate
        self.input_path = args.input_path
        self.output_path = args.output_path
        self.only_best_match = args.only_best_match
        self.database= args.database

    def validate_arguments(self):
        if self.mode == 'construct_database':
            for cur_path in self.input_path:
                if not os.path.exists(cur_path):
                    raise Exception('{} does not exist. `input_path` needs to contain fasta files that will be saved to DB'.format(cur_path))

def generate_tmp_file(size=10, chars=string.ascii_uppercase + string.digits): 
    filename = path.join('/tmp', ''.join(random.choice(chars) for x in range(size)))

    return filename



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

def get_combinations(filename, df_sequences, only_best_match=False, **kwargs):
    if only_best_match:
        return [get_best_matches(filename, df_sequences)]
    
    
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


## Creating fasta files that contains at most one gene for each species
def align_fasta_file(filename, **kwargs):
    try:
        #source = re.search(r"_(.*?)\.fasta", file).group(1)
        

        df_sequences = process_fasta_file(filename)
        if df_sequences.shape[0] <= 1:
            return False

        if len(df_sequences['species'].unique()) == df_sequences.shape[0]:
            combination = pd.Series(df_sequences.index, index=df_sequences['species'])
            combinations = list({combination.loc['Homo sapiens']: combination}.values())
        else:
            combinations = get_combinations(filename, df_sequences, **kwargs)        
        
        msa_results = {}

        for ind, combination in enumerate(combinations):
            fasta = align_combination(kwargs['output_dir'], combination, df_sequences)
            msa_results[ind] = fasta
        
        return df_sequences, combinations, msa_results
                        
    except Exception as e:
        print(e, filename)



if __name__ == '__main__':
    pipeline_options = Options()
    
    conn = connection_loader()
    sql_helper = MySqlHelper(conn)
    
    mode = pipeline_options.mode

    if mode == 'construct_database':
        input_dirs = pipeline_options.input_path
        output_dir = pipeline_options.output_path

        if pipeline_options.truncate == True:
            print("Truncating MSA tables")
            sql_helper.truncate_tables()
        else:
            print("Inserting the MSA tables without truncating them")

        #Directory which stores the processed fasta files which contains at most one gene for each species 
        processed_fasta_output_dir = path.join(output_dir, 'processed_fasta')
        
        if  os.path.exists(processed_fasta_output_dir):
            shutil.rmtree(processed_fasta_output_dir)
        
        os.makedirs(processed_fasta_output_dir)
        only_best_match = pipeline_options.only_best_match

        for input_dir in input_dirs:
            raw_fasta_files = [os.path.join(input_dir, filename) for filename in os.listdir(input_dir)]
            def chunk_callback(res, **kwargs):
                sql_helper.save_to_db(*res, only_best_match=only_best_match)

            results = run_parallel(
                                func = align_fasta_file, 
                                chunk_callback = chunk_callback,
                                args = raw_fasta_files, 
                                n_processes=26, 
                                only_best_match=only_best_match, 
                                output_dir=output_dir
                                )

    elif mode == 'single_fasta':
        filename = pipeline_options.input_path[0]
        if not os.path.exists(filename):
            print("File does not exist")

        print(create_alignment(filename))
