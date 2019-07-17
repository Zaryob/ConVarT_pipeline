import pandas as pd
import itertools
import os
from os import path
import sys
from pybiomart import Server
from tqdm import tqdm
from multiprocessing import Pool
import subprocess

species_to_dataset = pd.DataFrame([
	{'species_name':'Mus_musculus', 'table':'mapping_mouse', 'dataset':'mmusculus_gene_ensembl'},
	{'species_name':'Rattus_norvegicus', 'table':'mapping_rat', 'dataset':'rnorvegicus_gene_ensembl'},
	{'species_name':'Danio_rerio', 'table': 'mapping_zebrafish', 'dataset': 'drerio_gene_ensembl'},
	{'species_name':'Homo_sapiens', 'table':'mapping_human', 'dataset':'hsapiens_gene_ensembl'},
	{'species_name':'Pan_troglodytes', 'table': 'mapping_chimp', 'dataset':'ptroglodytes_gene_ensembl'},
	{'species_name':'Drosophila_melanogaster', 'table': 'mapping_fruitfly', 'dataset':'dmelanogaster_gene_ensembl'},
	{'species_name':'Macaca_mulatta', 'table': 'mapping_macaque', 'dataset': 'mmulatta_gene_ensembl'},
	{'species_name':'Xenopus_tropicalis', 'table': 'mapping_frog', 'dataset':'xtropicalis_gene_ensembl'},
	{'species_name':'Caenorhabditis_elegans', 'table': 'mapping_worm', 'dataset': 'celegans_gene_ensembl'}
])

biomart_data = pd.DataFrame()

server = Server(host='http://www.ensembl.org')

datasets = species_to_dataset['dataset'].tolist()
for dataset_name in datasets:
	dataset = (server.marts['ENSEMBL_MART_ENSEMBL']
                 .datasets[dataset_name])
	dataset_df = dataset.query(attributes=['ensembl_gene_id', 'entrezgene', 'description', 'external_gene_name'])

	biomart_data = pd.concat([biomart_data, dataset_df], axis=0)

allowed_db_list = ['Xenbase','Wormbase','HGNC','FLYBASE','EnsembleGenomes-Gn','MGI','ZFIN']

def process_gpff_list(f):
	gene_id_to_np =  {}
	db = {}
	other_dbs = set()
	for i, row in enumerate(f.readlines()):

		row = row.split(" ")
		np_id = row[0]
		data = row[1].replace("\n", "").split("=")
		data_key = data[0]
		data_value = [val.split(";") for val in 
						data[1].replace('"', '').replace(" ", "").split(",")]
		data_value = list(itertools.chain(*data_value))
		
		if np_id not in db:
			db[np_id] = {}
		
		if data_key == "gene":
			data_value = data_value[0]
		elif data_key == "locus_tag":
			data_value = [val.replace("CELE_", "") for val in data_value]
		
		elif data_key == "db_xrefs":
			data_value = {i.split(":")[0]:':'.join(i.split(":")[1:]) for i in data_value}
			other_dbs = other_dbs | set(data_value.keys())
			gene_id = data_value["GeneID"]

			if gene_id not in gene_id_to_np:
				gene_id_to_np[gene_id] = []

			gene_id_to_np[gene_id].append(np_id)

			data_value.pop("GeneID", None)
 		
		db[np_id][data_key] = data_value

	return db, gene_id_to_np, other_dbs

def imap_unordered_bar(func, args, n_processes = 10):
    p = Pool(n_processes)
    res_list = []
    with tqdm(total = len(args)) as pbar:
        for res in tqdm(p.imap_unordered(func, args)):
            pbar.update()
            res_list.append(res)
    pbar.close()
    p.close()
    p.join()
    return res_list

def gpff_to_list(gene_to_np):
	gene_id, np_ids = gene_to_np

	row = {}
	row['gene_id'] = [gene_id]
	row['protein_numbers'] = []
	row['gene_synonyms'] = []
	row['gene_symbol'] = []
	row['other_ids'] = []
	row['gene_description'] = []

	for np_id in np_ids:
		protein_data = db[np_id]

		if 'gene' in protein_data:
			row['gene_symbol'].append(protein_data['gene'])
			if 'locus_tag' in protein_data:
				row['gene_synonyms'] += protein_data['locus_tag']
		elif 'locus_tag' in protein_data:
			row['gene_symbol'] += protein_data['locus_tag']

		for key, value in protein_data['db_xrefs'].items():

			if key in allowed_db_list and value not in allowed_db_list:
				row['other_ids'].append(value)
			elif key == 'EnsemblGenomes-Tr':
				row['gene_synonyms'].append(value)

		if 'gene_synonym' in protein_data:
			row['gene_synonyms'] += protein_data['gene_synonym']

		row['protein_numbers'] += [np_id]
	
	biomart_row = biomart_data[biomart_data['NCBI gene ID']==int(gene_id)]
	biomart_row_by_symbol = biomart_data[biomart_data['Gene name']==row['gene_symbol'][0]]
	
	row['other_ids'] += biomart_row['Gene stable ID'].dropna().tolist()

	row['gene_description'] += biomart_row['Gene description'].dropna().tolist() + biomart_row_by_symbol['Gene description'].dropna().tolist()
	
	if len(row['gene_description']) > 1:
		row['gene_description'] = row['gene_description'][0:1]		

	row = {key: ','.join(list(set(value))) for key,value in row.items()}
	return row

gpff_dir = "/opt/current_project/db/mapping/gpff_work/gpff_list/"
csv_dir = "/opt/current_project/db/mapping/gpff_work/gpff_csv/"

df_all = pd.DataFrame()

if not path.exists(path.join(csv_dir, 'All.gpff.csv')):
	for species_file in os.listdir(gpff_dir):
		if '.list' not in species_file:
			continue


		with open(path.join(gpff_dir, species_file), 'r') as f:
			db, gene_id_to_np, other_dbs = process_gpff_list(f)

		rows = imap_unordered_bar(gpff_to_list, gene_id_to_np.items())		
		
		df = pd.DataFrame(rows)
		df['species_id'] = species_file.replace('.gpff.list', '').replace('_', ' ')
		df_all = pd.concat([df_all, df], axis=0)

	df_all.to_csv(path.join(csv_dir, 'All.gpff.csv'), sep=';', index=False)
else:
	df_all = pd.read_csv(path.join(csv_dir, 'All.gpff.csv'), sep=';')
import pymysql 
con = pymysql.connect(host='127.0.0.1', unix_socket='/opt/lampp/var/mysql/mysql.sock', 
            user='root', passwd='', db='current_project')
cur = con.cursor()

for index, row in df_all.iterrows():
	gene_id = row['gene_id']

	for key, values in row.to_dict().items():
		if key == 'gene_id':
			continue
		if key[-1] == 's':
			key = key[:-1]
		try:
			if values == 'nan':
				continue

			if key == 'gene_description':
				cur.execute("INSERT INTO ncbi_gene_meta(ncbi_gene_id, meta_key, meta_value) VALUES(%s, %s, %s)", (gene_id, key, values))
				continue

			for value in str(values).split(','):
				if value == 'nan':
					continue
				cur.execute("INSERT INTO ncbi_gene_meta(ncbi_gene_id, meta_key, meta_value) VALUES(%s, %s, %s)", (gene_id, key, value))
		except Exception as e:
			print(e, gene_id, gene_id, key, value)

	con.commit()

#subprocess.run("mysql -u root current_project < /opt/current_project/pipeline/sql/insert_gpff_files.sql", shell=True)