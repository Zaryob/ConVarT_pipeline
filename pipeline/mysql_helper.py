import re 
import hashlib

class MySqlHelper():

    def __init__(self, conn):

        self.cur, self.con = conn.cursor(), conn

    def write_msa_to_db(self, fasta, combination):

        self.cur.execute("INSERT INTO msa (fasta, alignment_method) VALUES(%s, %s)", (fasta, 'clustalw'))
        
        msa_id = self.cur.lastrowid

        for np_id in combination.values:
            np_id_without_version = np_id.split('.')[0]

            self.cur.execute("SELECT convart_gene_id FROM convart_gene_to_db WHERE db_id = %s LIMIT 1", 
                                            (np_id_without_version))
            row = self.cur.fetchone()
            convart_gene_id = row[0]
            self.con.commit()

            self.cur.execute("INSERT INTO msa_gene (msa_id, convart_gene_id) VALUES(%s, %s)", (msa_id, convart_gene_id))
        
        self.con.commit()

        return msa_id

    def write_fasta_to_db(self, fasta, db_id, species):
        
        db_id_version = None
        if len(db_id.split(".")) == 2:
            db_id, db_id_version = db_id.split('.')[0], db_id.split('.')[1]
        elif len(db_id.split(".")) > 2:
            db_id, db_id_version = ".".join(db_id.split(".")[0:2]), db_id.split('.')[-1]

        if db_id.startswith('NP'):
            db = 'NP'
        elif db_id.startswith('ENST'):
            db = 'ENST'
        elif db_id.startswith('ENSMUST'):
            db = 'ENSMUST'
        else:
            db = 'OTHER'

        fasta = re.sub(r'>(.*?)\n', '', fasta).replace('\n','').replace("'", '')
        seq_hash = hashlib.md5((fasta+species).encode('utf-8')).hexdigest()

        self.cur.execute("SELECT id FROM convart_gene WHERE hash=%s and species_id=%s", (seq_hash, species))
        
        if self.cur.rowcount == 0:
            self.cur.execute("INSERT IGNORE INTO convart_gene (sequence, species_id, hash) VALUES(%s, %s, %s)", 
                                (fasta, species, seq_hash))
            
            convart_gene_id = self.con.insert_id()
            self.con.commit()
        else:
            row = self.cur.fetchone()
            convart_gene_id = row[0]
        self.cur.execute("INSERT IGNORE INTO convart_gene_to_db (convart_gene_id, db, db_id, db_id_version) VALUES(%s, %s, %s, %s)", 
                        (convart_gene_id, db, db_id, db_id_version))

        return convart_gene_id

    def write_best_combination_to_db(self, msa_id, convart_gene_id):

        self.cur.execute("INSERT IGNORE INTO msa_best_combination (msa_id, convart_gene_id) VALUES(%s, %s)", 
                                (msa_id, convart_gene_id))
        self.con.commit()


    def save_to_db(self, df_sequences, combinations, msa_results, only_best_match=False):
        for gene_id, row in df_sequences.iterrows():
            
            convart_gene_id = self.write_fasta_to_db(row['sequence'], gene_id, row['species'])
            if row['species'] == 'Homo sapiens':
                human_convart_gene_id = convart_gene_id
            
        for ind, combination in enumerate(combinations):
            fasta = msa_results[ind]
            msa_id = self.write_msa_to_db(fasta, combination)
            if only_best_match:
                self.write_best_combination_to_db(msa_id, human_convart_gene_id)

    def truncate_tables(self):
        cur, con = self.cur, self.con
        
        cur.execute('TRUNCATE TABLE msa')
        cur.execute('TRUNCATE TABLE msa_gene')
        cur.execute('TRUNCATE TABLE msa_best_combination')
        cur.execute('TRUNCATE TABLE convart_gene')
        cur.execute('TRUNCATE TABLE convart_gene_to_db')
        con.commit()