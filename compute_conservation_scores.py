import sys
from sqlalchemy import create_engine
import pandas as pd
import os
from os import path
import re
import pymysql
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import warnings
warnings.filterwarnings("ignore")
from tqdm import tqdm
from multiprocessing import Pool
def imap_unordered_bar(func, args, number_of_elements, n_processes = 30):
    scores = []
    p = Pool(n_processes)

    with tqdm(total = number_of_elements) as pbar:
        for res in tqdm(p.imap_unordered(func, args, chunksize=n_processes)):
            scores += res
            pbar.update()

    pbar.close()
    p.close()
    p.join()

    return scores

engine = create_engine("mysql+pymysql://{user}:{pw}@localhost/{db}"
                       .format(user="root",
                               pw="",
                               db="current_project"))


plot_dir = '/opt/current_project/tmp/clinvar/'
alignment_dir = "/opt/current_project/results/alignments/"

polarity = {
	'A': 'nonpolar',
	'R': 'polar',
	'N': 'polar',
	'D': 'polar',
	'C': 'nonpolar',
	'E': 'polar',
	'Q': 'polar',
	'G': 'nonpolar',
	'H': 'polar',
	'I': 'nonpolar',
	'L': 'nonpolar',
	'K': 'polar',
	'M': 'nonpolar',
	'F': 'nonpolar',
	'P': 'nonpolar',
	'S': 'polar',
	'T': 'polar',
	'W': 'nonpolar',
	'Y': 'polar',
	'V': 'nonpolar',
	'U': 'polar',
	'-': '-'
}
def sqlize(significance_list):
	return "'"+ "', '".join(significance_list) + "'"

con = pymysql.connect(host='127.0.0.1', unix_socket='/opt/lampp/var/mysql/mysql.sock', user='root', passwd='', db='current_project')
cur = con.cursor()

ncbi_to_transcript_id = pd.read_sql("SELECT gdb.db_id AS NCBI, gdb.convart_gene_id AS id, ncbi.meta_value AS gene_symbol "
								"FROM convart_gene_to_db AS gdb " +
								"LEFT JOIN ncbi_gene_meta AS ncbi ON ncbi.ncbi_gene_id = gdb.db_id WHERE gdb.db='NCBI' "+
								"AND ncbi.meta_key = 'gene_symbol' GROUP BY gdb.convart_gene_id", con)
ncbi_to_transcript_id.set_index('id', drop=False, inplace=True)


genes = pd.read_sql("SELECT DISTINCT g.id, gdb_enst.db_id AS ENST, gdb_nm.db_id AS NM, m.fasta, m.id as msa_id, g.sequence "+
			"FROM convart_gene AS g "+
			"LEFT JOIN convart_gene_to_db AS gdb_enst ON gdb_enst.convart_gene_id=g.id AND gdb_enst.db='ENST' " + 
			"LEFT JOIN convart_gene_to_db AS gdb_nm ON gdb_nm.convart_gene_id=g.id AND gdb_nm.db='NM' " + 
			"INNER JOIN msa_best_combination AS mb ON mb.convart_gene_id = g.id " +
			"INNER JOIN msa AS m ON m.id = mb.msa_id "
			"WHERE g.species_id='Homo sapiens'", con)


enst_genes = genes['id'].str.startswith('ENST')
genes.loc[enst_genes, 'ENST'] = genes[enst_genes]['id'].str.split('.').apply(lambda x: x[0])

benign_list = sqlize(["Benign","Benign, association","Benign/Likely benign","Benign/Likely benign, Affects","Benign/Likely benign, association",
			"Benign/Likely benign, drug response","Benign/Likely benign, drug response, risk factor","Benign/Likely benign, other","Benign/Likely benign, protective",
			"Benign/Likely benign, protective, risk factor","Benign/Likely benign, risk factor","Benign, other","Benign, risk factor","Likely benign",
			"Likely benign, drug response, other","Likely benign, other","Likely benign, risk factor"])
vus_list = sqlize(["Uncertain significance","Uncertain significance, drug response","Uncertain significance, other","Uncertain significance, risk factor"])
pathogenic_list = sqlize(["Pathogenic","Pathogenic, Affects","Pathogenic, association, protective","Pathogenic, drug response","Pathogenic/Likely pathogenic",
"Pathogenic/Likely pathogenic, drug response","Pathogenic/Likely pathogenic, other","Pathogenic/Likely pathogenic, risk factor","Pathogenic, other",
"Pathogenic, other, risk factor","Pathogenic, protective","Pathogenic, risk factor","Likely pathogenic","Likely pathogenic, association",
"Likely pathogenic, drug response","Likely pathogenic, other","Likely pathogenic, risk factor"])


clinvar_data = pd.read_sql("SELECT gdb.convart_gene_id AS transcript_id, GROUP_CONCAT(DISTINCT position) as positions, COUNT(DISTINCT position) "
			+"as variation_count FROM clinvar INNER JOIN convart_gene_to_db AS gdb ON gdb.db_id=nm_id WHERE position > 0 "+
			" GROUP BY gdb.convart_gene_id", con)
pathogenic_clinvar_data = pd.read_sql("SELECT gdb.convart_gene_id AS transcript_id, GROUP_CONCAT(DISTINCT position) as positions, COUNT(DISTINCT position) "
			+"as variation_count FROM clinvar INNER JOIN convart_gene_to_db AS gdb ON gdb.db_id=nm_id WHERE position > 0 AND clinical_significance IN ("+pathogenic_list+")"+
			" GROUP BY gdb.convart_gene_id", con)
vus_clinvar_data = pd.read_sql("SELECT gdb.convart_gene_id AS transcript_id, GROUP_CONCAT(DISTINCT position) as positions, COUNT(DISTINCT position) "
			+"as variation_count FROM clinvar INNER JOIN convart_gene_to_db AS gdb ON gdb.db_id=nm_id WHERE position > 0 AND clinical_significance IN ("+vus_list+")"+
			" GROUP BY gdb.convart_gene_id", con)
benign_clinvar_data = pd.read_sql("SELECT gdb.convart_gene_id AS transcript_id, GROUP_CONCAT(DISTINCT position) as positions, COUNT(DISTINCT position) "
			+"as variation_count FROM clinvar INNER JOIN convart_gene_to_db AS gdb ON gdb.db_id=nm_id WHERE position > 0 AND clinical_significance IN ("+benign_list+")"+
			" GROUP BY gdb.convart_gene_id", con)
gnomad_data = pd.read_sql("SELECT gnomad.canonical_transcript AS transcript_id, GROUP_CONCAT(DISTINCT position) AS positions, COUNT(DISTINCT position) "
					+"AS variation_count FROM gnomad WHERE position > 0 GROUP BY gnomad.canonical_transcript", con)

combined_clinvar_gnomad_data = []


variations = {
	'variant': pd.concat([clinvar_data, gnomad_data], axis=0, ignore_index=True),
	'pathogenic': pathogenic_clinvar_data,
	'vus': vus_clinvar_data,
	'benign': benign_clinvar_data
}


def conservation_scores_by_gene(row):
	scores = []

	ind, msa_details = row

	msa_details = msa_details.to_dict()
	sequences = ('\n' + msa_details['fasta']).split('\n>')

	sequences = [{
			'transcript_id': re.sub('^(.*?) (.*?)$', r'\1', seq.split('\n')[0]),
			'species_id': re.sub(r'^(.*?)\[([a-zA-Z0-9 ]+?)\]$', r'\2', seq.split('\n')[0]),
			'fasta': re.sub('^(.*?)\n', '', seq).replace('\n','').replace(' ', '')} for seq in sequences]
	sequences = list(filter(lambda x : x['fasta'].replace('\n','') !='', sequences))

	human_seq = sequences[0]['fasta']
	seq_length = len(human_seq)
	n_sequences = len(sequences)
	
	variation_positions_by_data = {}

	for key, data in variations.items():
		variation_positions_by_data[key] = set()

		if data[data['transcript_id'] == msa_details['id']].shape[0] > 0:
			cur_positions = data[data['transcript_id'] == msa_details['id']]['positions'].values[0].split(',')
			[variation_positions_by_data[key].add(int(cur_position)) for cur_position in cur_positions if cur_position != '']
			
	for pairwise_id in range(1, len(sequences)):
		exact_match_score = 0
		polarity_match_score = 0
		
		exact_scores_by_data = {}
		polarity_scores_by_data = {}

		for key, data in variations.items():
			exact_scores_by_data[key] = 0
			polarity_scores_by_data[key] = 0

		pair = sequences[pairwise_id]
		pair_seq = pair['fasta']
		human_seq_position = 0
		
		for i in range(seq_length):
			if len(human_seq) != len(pair_seq):
				print(len(human_seq), len(pair_seq), msa_details['msa_id'])
				sys.exit(-1)
			if human_seq[i] == '-':
				continue
			human_seq_position += 1

			if human_seq[i] == 'X' or pair_seq[i] == 'X':
				continue

			aligned_aminoacids = 1 if len({seq[i] for seq in [human_seq, pair_seq]})==1 else 0
			exact_match_score += aligned_aminoacids

			aligned_aminoacids_polarity = 1 if  len({polarity[seq[i]] for seq in [human_seq, pair_seq]}) == 1 else 0
			polarity_match_score += aligned_aminoacids_polarity

			for key, positions in variation_positions_by_data.items():
				if human_seq_position not in positions:
					continue
				exact_scores_by_data[key] += aligned_aminoacids
				polarity_scores_by_data[key] += aligned_aminoacids_polarity

		scores.append({'transcript_id': msa_details['id'], 'pairwise_transcript_id': pair['transcript_id'],
						'pairwise_species_id': pair['species_id'],
						'score_type': 'exact_match_score', 'score': exact_match_score, 
						'number_of_aminoacids': len(msa_details['sequence'])})
		scores.append({'transcript_id': msa_details['id'], 'pairwise_transcript_id': pair['transcript_id'],
						'pairwise_species_id': pair['species_id'],
						'score_type': 'polarity_match_score', 'score': polarity_match_score,
						'number_of_aminoacids': len(msa_details['sequence'])})

		#Variant scores
		for key, exact_score in exact_scores_by_data.items():
			scores.append({'transcript_id': msa_details['id'], 'pairwise_transcript_id': pair['transcript_id'],
						'pairwise_species_id': pair['species_id'],
						'score_type': key+'_exact_match_score', 'score': exact_score, 
						'number_of_aminoacids': len(variation_positions_by_data[key])})
		for key, polarity_score in polarity_scores_by_data.items():
			scores.append({'transcript_id': msa_details['id'], 'pairwise_transcript_id': pair['transcript_id'],
						'pairwise_species_id': pair['species_id'],
						'score_type': key+'_polarity_match_score', 'score': polarity_score, 
						'number_of_aminoacids': len(variation_positions_by_data[key])})

	return scores

scores = imap_unordered_bar(conservation_scores_by_gene, genes.iterrows(), genes.shape[0], 30)

df = pd.DataFrame(scores)

def get_id(x):
	global ncbi_to_transcript_id
	search_params = [x]
	if x.startswith('ENST'):
		search_params += [x.split('.')[0]]

	for param in search_params:
		if param in ncbi_to_transcript_id.index:
			return ncbi_to_transcript_id.loc[param, 'gene_symbol']

	return None

df.loc[:, 'gene_symbol'] = df['transcript_id'].apply( get_id )

column_map = {
	'number_of_aminoacids': 'aminoacid_number',
	'pairwise_species_id': 'specie',
	'pairwise_transcript_id':'transcript_id_specie',
	'score': 'score',
	'score_type': 'score_type',
	'transcript_id': 'transcript_id',
	'gene_symbol': 'gene_symbol'
}

df.columns = [column_map[col] for col in df.columns]
df.loc[df['gene_symbol'].isna(), 'gene_symbol'] = ''
cur.execute('TRUNCATE TABLE conservation_scores')
con.commit()
engine = create_engine("mysql+pymysql://{user}:{pw}@localhost/{db}"
                       .format(user="root",
                               pw="",
                               db="current_project"))

df.to_sql('conservation_scores', engine, if_exists='append', index=False)

#df.to_csv('/opt/current_project/results/statistics/conservation_scores.csv', index=None)
