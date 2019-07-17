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

def parallel_loop(func, args, n_processes = 30):
    p = Pool(n_processes)

    with tqdm(total = len(args)) as pbar:
        for res in tqdm(p.imap_unordered(func, args, chunksize=n_processes)):
            pbar.update()

    pbar.close()
    p.close()
    p.join()

with open('../config.yml') as f:
    config = yaml.load(f)
project_path = config['PROJECT']
db_path = config['DATABASE']

nm_to_transcript_id = pd.read_csv(path.join(db_path, 'mapping', 'NM_NP_GeneID.list'), sep='\t', header=None)
nm_to_transcript_id.columns = ['retrieval_id', 'transcript_id']
enst_to_transcript_id = pd.read_csv(path.join(db_path, 'mapping', 'ENST_ENSTV_GeneID.list'), header=None)
enst_to_transcript_id.columns = ['retrieval_id', 'transcript_id']

retrieval_id_to_transcript_id = pd.concat([nm_to_transcript_id, enst_to_transcript_id], axis=0, ignore_index=True)
retrieval_id_to_transcript_id = retrieval_id_to_transcript_id.set_index('retrieval_id')


## Homology parse

human_transcripts = glob(path.join(project_path, 'results', 'seqs*', '*'))
homology_list = pd.read_csv(path.join(db_path,'orthology','Homology.list'), sep='\t', dtype=str)

species_list = ['Pan_troglodytes', 'Macaca_mulatta', 'Rattus_norvegicus', 'Mus_musculus', 'Danio_rerio', 
            'Xenopus_tropicalis', 'Drosophila_melanogaster', 'Caenorhabditis_elegans']
homology_list.columns = ['Homo_sapiens'] + species_list
homology_list = homology_list.set_index('Homo_sapiens')

transcript_lists = {}
fasta_lists = {}

enst_to_gene_id = pd.read_csv(path.join(db_path, 'mapping', 'NewCurated_ENSTvsGENEID.csv'))
enst_to_gene_id.columns = ['transcript_id', 'gene_id']
np_to_gene_id = pd.read_csv(path.join(db_path, 'mapping', 'NM_NP_GeneID.list'), sep='\t').iloc[:, 1:] 
np_to_gene_id.columns = ['transcript_id', 'gene_id']

transcript_id_to_gene_id = pd.concat([enst_to_gene_id, np_to_gene_id], axis=0, ignore_index=True)

gene_id_to_transcript_ids = pd.concat([enst_to_gene_id, np_to_gene_id], axis=0, ignore_index=True)

for species in species_list:
    species_list_file = path.join(db_path, 'mapping', 'NP_tables', species+'.list')
    print(species_list_file)
    transcript_lists[species] = pd.read_csv(species_list_file, sep='\t')
    transcript_list = transcript_lists[species][['Protein product', 'GeneID']]
    transcript_list.columns = ['transcript_id', 'gene_id']
    
    transcript_id_to_gene_id = pd.concat([transcript_id_to_gene_id, transcript_list],
                                        axis=0, ignore_index=True)
    gene_id_to_transcript_ids = pd.concat([gene_id_to_transcript_ids, transcript_list],
                                        axis=0, ignore_index=True)
    
    fasta_lists[species] = {} 
    fasta_sequence_file_gz = path.join(db_path, 'proteins', species+'.faa.gz')
    fasta_sequence_file = path.join(db_path, 'proteins', species+'.faa')

    with gzip.open(fasta_sequence_file_gz, 'rb') as handle:
        fasta = handle.read()
        with open(fasta_sequence_file, 'wb') as f_in:
            f_in.write(fasta)
    
    fasta_sequences = SeqIO.parse(fasta_sequence_file, "fasta")
    for fasta in tqdm(fasta_sequences):
        transcript_id, sequence = fasta.id, str(fasta.seq)
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


def generate_fasta(transcript_file):
    with open(transcript_file, 'r') as f:
        current_fasta = f.read()
    source = '_' + transcript_file.split('_')[-1]
    human_transcript_id_orig = current_fasta.split('\n')[0].split(' ')[0].replace('>', '') 
    human_transcript_id = human_transcript_id_orig.split('.')[0]
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
                
                current_fasta += '>'+transcript_id+' ['+species.replace('_',' ')+']'    
                current_fasta += '\n'+fasta_lists[species][transcript_id]
                current_fasta += '\n'
                more_than_one = True

    if not more_than_one:
        return False
    with open(path.join(project_path+'/results/trial_seqs/'+\
                        human_transcript_id_orig + source), 'w') as f:
        f.write(current_fasta)

if path.exists(project_path+'/results/trial_seqs/'):
    shutil.rmtree(project_path+'/results/trial_seqs/')
os.makedirs(project_path+'/results/trial_seqs/')

parallel_loop(generate_fasta, human_transcripts)
