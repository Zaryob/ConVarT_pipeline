import pandas as pd
from os import path
from Bio import SeqIO
import gzip
from tqdm import tqdm
from glob import glob
from multiprocessing import Pool
import shutil
import os
import yaml
import argparse

SPECIES_LIST = ['Pan_troglodytes', 'Macaca_mulatta', 'Rattus_norvegicus', 'Mus_musculus', 'Danio_rerio', 
            'Xenopus_tropicalis', 'Drosophila_melanogaster', 'Caenorhabditis_elegans']

def create_formatted_fasta(transcript_id, species, sequence):
    current_fasta = '>'+transcript_id+' ['+species.replace('_',' ')+']'    
    current_fasta += '\n'+sequence
    current_fasta += '\n'

    return current_fasta

def parallel_loop(func, args, n_processes = 30):
    p = Pool(n_processes)

    with tqdm(total = len(args)) as pbar:
        for res in tqdm(p.imap_unordered(func, args, chunksize=n_processes)):
            pbar.update()

    pbar.close()
    p.close()
    p.join()

class Options():
    
    def __init__(self):
       self.parse_arguments()      
       self.validate_arguments()
       self.parse_config()
        
    def parse_config(self):
        dir_path = os.path.dirname(os.path.realpath(__file__))

        with open(dir_path+'/../config.yml') as f:
            config = yaml.load(f)

        self.project_path = config['PROJECT']
        
        self.db_path = config['DATABASE']
        self.output_path = self.project_path+'/results/{}/'.format(self.output_path)

        
    def parse_arguments(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("--included_species", type=str, nargs='+', required=True,
                            help="List of species that are included in the fasta files that will be created")
        parser.add_argument("--output_path", type=str, default='trial_seqs', 
                        help="the subdirectory that outputs will be saved in $PROJECT_PATH/results/{output_path}, where PROJECT_PATH is set in config.yml. ")
        args = parser.parse_args()
        self.args = args
        self.included_species = args.included_species
        self.output_path = args.output_path
    
    def validate_arguments(self):
        included_species = []
        for species in self.included_species:
            processed_species = species.capitalize().replace(' ', '_')
            
            if processed_species not in SPECIES_LIST:
                raise Exception("{} not found in the available ortholog species list".format(species))
            else:
                included_species.append(processed_species)
        
        self.included_species = included_species

if __name__ == '__main__':
    pipeline_options = Options()

    nm_to_transcript_id = pd.read_csv(path.join(pipeline_options.db_path, 'mapping', 'NM_NP_GeneID.list'), sep='\t', header=None).iloc[:, :2]
    nm_to_transcript_id.columns = ['retrieval_id', 'transcript_id']

    ## Homology parse

    homology_list = pd.read_csv(path.join(pipeline_options.db_path,'orthology','Homology.list'), sep='\t', dtype=str)


    homology_list.columns = ['Homo_sapiens'] + SPECIES_LIST
    homology_list = homology_list.set_index('Homo_sapiens')

    transcript_lists = {}
    fasta_lists = {}

    enst_to_gene_id = pd.read_csv(path.join(pipeline_options.db_path, 'mapping', 'NewCurated_ENSTvsGENEID.csv'))
    enst_to_gene_id.columns = ['transcript_id', 'gene_id']
    np_to_gene_id = pd.read_csv(path.join(pipeline_options.db_path, 'mapping', 'NM_NP_GeneID.list'), sep='\t').iloc[:, 1:] 
    np_to_gene_id.columns = ['transcript_id', 'gene_id']

    transcript_id_to_gene_id = pd.concat([enst_to_gene_id, np_to_gene_id], axis=0, ignore_index=True)

    gene_id_to_transcript_ids = pd.concat([enst_to_gene_id, np_to_gene_id], axis=0, ignore_index=True)

    for species in pipeline_options.included_species + ['Homo_sapiens']:
        species_list_file = path.join(pipeline_options.db_path, 'mapping', 'NP_tables', species+'.list')

        transcript_lists[species] = pd.read_csv(species_list_file, sep='\t')
        transcript_list = transcript_lists[species][['Protein product', 'GeneID']]
        transcript_list.columns = ['transcript_id', 'gene_id']
        if species != 'Homo_sapiens':
            transcript_id_to_gene_id = pd.concat([transcript_id_to_gene_id, transcript_list],
                                                axis=0, ignore_index=True)
            gene_id_to_transcript_ids = pd.concat([gene_id_to_transcript_ids, transcript_list],
                                                axis=0, ignore_index=True)
        
        fasta_lists[species] = {} 
        fasta_sequence_file_gz = path.join(pipeline_options.db_path, 'proteins', species+'.faa.gz')
        fasta_sequence_file = path.join(pipeline_options.db_path, 'proteins', species+'.faa')

        if not os.path.exists(fasta_sequence_file):
            with gzip.open(fasta_sequence_file_gz, 'rb') as handle:
                fasta = handle.read()
                with open(fasta_sequence_file, 'wb') as f_in:
                    f_in.write(fasta)
        
        fasta_sequences = SeqIO.parse(fasta_sequence_file, "fasta")
        for fasta in tqdm(fasta_sequences):
            transcript_id, sequence = fasta.id, str(fasta.seq)

            if species != 'Homo_sapiens':
                transcript_id = transcript_id.split('.')[0]
    
            fasta_lists[species][transcript_id] = sequence

    transcript_id_to_gene_id = transcript_id_to_gene_id.dropna()
    transcript_id_to_gene_id.loc[:, 'transcript_id'] = transcript_id_to_gene_id['transcript_id'].apply(lambda x: x.split('.')[0])
    transcript_id_to_gene_id.loc[:, 'gene_id'] = transcript_id_to_gene_id.loc[:, 'gene_id'].astype(str) 
    transcript_id_to_gene_id = transcript_id_to_gene_id.set_index('transcript_id')

    gene_id_to_transcript_ids = gene_id_to_transcript_ids.dropna()
    gene_id_to_transcript_ids.loc[:, 'gene_id'] = gene_id_to_transcript_ids.loc[:, 'gene_id'].astype(str)

    gene_id_to_transcript_ids.loc[:, 'transcript_id'] = gene_id_to_transcript_ids['transcript_id'].apply(lambda x: x.split('.')[0])
    gene_id_to_transcript_ids = gene_id_to_transcript_ids.groupby('gene_id')['transcript_id'].apply(lambda x: ','.join(x))


    def generate_fasta(human_transcript_id_with_version):
        
        current_fasta = create_formatted_fasta(human_transcript_id_with_version, 'Homo_sapiens', fasta_lists['Homo_sapiens'][human_transcript_id_with_version])
        human_transcript_id = human_transcript_id_with_version.split('.')[0]

        if human_transcript_id not in transcript_id_to_gene_id.index:
            return False
            
        if pd.isna(transcript_id_to_gene_id.loc[human_transcript_id, 'gene_id']):
            return False
        if ' ' in str(transcript_id_to_gene_id.loc[human_transcript_id, 'gene_id']):
            print(human_transcript_id, 
                transcript_id_to_gene_id.loc[human_transcript_id, 'gene_id'])
            return False
        
        human_gene_id = transcript_id_to_gene_id.loc[human_transcript_id, 'gene_id']
        
        if human_gene_id not in homology_list.index:
            return False
        
        saved_ids = set()
        more_than_one = False
        for species, homolog_gene_ids in homology_list.loc[human_gene_id].to_dict().items():
            if species not in pipeline_options.included_species:
                continue
            
            if pd.isna(homolog_gene_ids):
                continue
            homolog_gene_ids = homolog_gene_ids.split(',')
            for gene_id in homolog_gene_ids:
                if gene_id == '':
                    continue
                if gene_id not in gene_id_to_transcript_ids.index:
                    continue

                for transcript_id in gene_id_to_transcript_ids.loc[gene_id].split(','):
                    if transcript_id not in fasta_lists[species]:
                        continue
                    if transcript_id in saved_ids:
                        continue
                    saved_ids.add(transcript_id)
                    
                    current_fasta += create_formatted_fasta(transcript_id, species, sequence)
                    
                    more_than_one = True

        if not more_than_one:
            return False

        with open(path.join(pipeline_options.output_path+human_transcript_id_with_version), 'w') as f:
            f.write(current_fasta)


    if path.exists(pipeline_options.output_path):
        shutil.rmtree(pipeline_options.output_path)
    
    os.makedirs(pipeline_options.output_path)
    
    parallel_loop(generate_fasta, fasta_lists['Homo_sapiens'].keys())
