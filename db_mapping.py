import sys
import pymysql
import pandas as pd
import warnings
from tqdm import tqdm
warnings.simplefilter(action="ignore", category=pymysql.Warning)

def save_mapping_to_db(df_basics, con, cur):
	for np_id, row in df_basics.iterrows():
		for db, db_id in row.items():
			cur.execute("INSERT IGNORE INTO convart_gene_to_db (convart_gene_id, db, db_id) VALUES (%s, %s, %s)", 
							(np_id, db, db_id))

if __name__ == '__main__':
    con = pymysql.connect(host='127.0.0.1', unix_socket='/opt/lampp/var/mysql/mysql.sock', 
            user='root', passwd='', db='current_project')
    cur = con.cursor()

    df_basics = pd.read_csv('/opt/current_project/db/mapping/NM_NP_GeneID.list', header=None, sep='\t')
    df_basics.columns = ['NM', 'NP', 'NCBI']
    for ind, row in tqdm(df_basics.iterrows()):
        cur.execute("UPDATE clinvar SET np_id=%s WHERE nm_id=%s", 
            (row['NP'].split('.')[0], row['NM'].split('.')[0]))
    con.commit()
               
    #save_mapping_to_db(df_basics, con, cur)
    sys.exit(0)
    
    #df_duplicates = pd.read_sql("SELECT msa.id FROM msa INNER JOIN msa_gene AS mg ON mg.msa_id=msa.id WHERE SUBSTRING_INDEX(mg.convart_gene_id, '.', 1) "+
    #				"IN (SELECT db_id FROM convart_gene_to_db WHERE db = 'ENST')", con)
    #duplicate_ensts = [str(i) for i in df_duplicates.iloc[:, 0].tolist()]
    #sqlized_duplicate_ensts = '"' + '", "'.join(duplicate_ensts)+'"'

    #cur.execute('DELETE FROM msa WHERE id IN ('+sqlized_duplicate_ensts+')')
    #cur.execute('DELETE FROM msa_best_combination WHERE msa_id IN ('+sqlized_duplicate_ensts+')')
    #cur.execute('DELETE FROM msa_gene WHERE msa_id IN ('+sqlized_duplicate_ensts+')')
    #con.commit()



    con.commit()

    df_basics = pd.read_csv('/opt/current_project/db/mapping/MissingGenes.csv', sep=',')
    df_basics = df_basics[['No Matching', 'NCBI', 'Gene Name']]
    df_basics = df_basics.set_index('No Matching')
    df_basics = df_basics.dropna().drop_duplicates()

    for transcript_id, row in df_basics.iterrows():
        if transcript_id.startswith('ENST'):
            transcript_id = transcript_id.split('.')[0]

        cur.execute("INSERT IGNORE INTO ncbi_gene_meta (ncbi_gene_id, meta_key, meta_value) VALUES (%s, 'protein_number', %s)", 
                            (row['NCBI'], transcript_id))

        cur.execute("INSERT IGNORE INTO ncbi_gene_meta (ncbi_gene_id, meta_key, meta_value) VALUES (%s, 'gene_symbol', %s)", 
                            (row['NCBI'], row['Gene Name']))
    con.commit()
    df_basics.drop('Gene Name', axis=1, inplace=True)
    save_mapping_to_db(df_basics, con, cur)


    df_basics = pd.read_csv('/opt/current_project/db/mapping/NewCurated_ENSTvsGENEID.csv', sep=',')
    df_basics.columns = ['ENST', 'NCBI']
    df_basics = df_basics.set_index('ENST')
    df_basics = df_basics.dropna().drop_duplicates()

    for np_id, row in df_basics.iterrows():

        cur.execute("INSERT IGNORE INTO ncbi_gene_meta (ncbi_gene_id, meta_key, meta_value) VALUES (%s, 'protein_number', %s)", 
                            (row['NCBI'], np_id))
    con.commit()

    for np_id, row in df_basics.iterrows():

        cur.execute("INSERT IGNORE INTO ncbi_gene_meta (ncbi_gene_id, meta_key, meta_value) VALUES (%s, 'species_id', %s)", 
                            (row['NCBI'], 'Homo sapiens'))


    con.commit()

    save_mapping_to_db(df_basics, con, cur)


    df_basics = pd.read_csv('/opt/current_project/db/mapping/SequenceMatch_ClinVar_gnomAD.list', sep=',')
    df_basics = df_basics[df_basics['needle_score']=='100.0%']
    df_basics = df_basics[['seq_from_clinvar', 'seq_from_gnomAD']]
    df_basics.columns = ['NP', 'ENST']
    df_basics = df_basics.set_index('NP')

    save_mapping_to_db(df_basics, con, cur)


