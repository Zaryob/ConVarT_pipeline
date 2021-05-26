import pandas as pd
import itertools
import os
import sys
from tqdm import tqdm
from multiprocessing import Pool
import urllib.request
import pathlib
from os import path
import glob
from pipeline.parallel_helper import run_parallel
from pipeline.commons import connection_loader, config_loader, catched_subprocess, BiomartHelper


class GpffHelper:
    gpff_links = [
            {
                "species_name": "Homo_sapiens",
                "link": "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_protein.gpff.gz",
            },
            {
                "species_name": "Pan_troglodytes",
                "link": "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/880/755/GCF_002880755.1_Clint_PTRv2/GCF_002880755.1_Clint_PTRv2_protein.gpff.gz",
            },
            {
                "species_name": "Macaca_mulatta",
                "link": "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/772/875/GCF_000772875.2_Mmul_8.0.1/GCF_000772875.2_Mmul_8.0.1_protein.gpff.gz",
            },
            {
                "species_name": "Rattus_norvegicus",
                "link": "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/895/GCF_000001895.5_Rnor_6.0/GCF_000001895.5_Rnor_6.0_protein.gpff.gz",
            },
            {
                "species_name": "Mus_musculus",
                "link": "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_protein.gpff.gz",
            },
            {
                "species_name": "Danio_rerio",
                "link": "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/035/GCF_000002035.6_GRCz11/GCF_000002035.6_GRCz11_protein.gpff.gz",
            },
            {
                "species_name": "Xenopus_tropicalis",
                "link": "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/195/GCF_000004195.3_Xenopus_tropicalis_v9.1/GCF_000004195.3_Xenopus_tropicalis_v9.1_protein.gpff.gz",
            },
            {
                "species_name": "Drosophila_melanogaster",
                "link": "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_protein.gpff.gz",
            },
            {
                "species_name": "Caenorhabditis_elegans",
                "link": "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/GCF_000002985.6_WBcel235/GCF_000002985.6_WBcel235_protein.gpff.gz",
            },
        ]
    

    @staticmethod
    def download_files(download_path):
        for gpff in GpffHelper.gpff_links:
            species, link = gpff["species_name"], gpff["link"]
            filename = f"{download_path}/{species}.gpff.gz"
            if not os.path.exists(filename):
                urllib.request.urlretrieve(link, filename)

    @staticmethod
    def gpff_to_list(gpff_path, output_path):

        dir_path = os.path.dirname(os.path.realpath(__file__))

        for gpff_file in glob.glob(path.join(gpff_path, "*.gz")):
    
            output_file = path.join(output_path, gpff_file.split("/")[-1].replace('.gz', '') +'.list')
            catched_subprocess(
                f"""bash -c "awk -f {dir_path}/gpff_ret.awk <(zcat {gpff_file}) > {output_file}" """
            )
    
    @staticmethod
    def process_gpff_list(gpff_list_file):
        gpff_list = []
        with open(gpff_list_file, 'r') as f:
            gpff_list = f.readlines()

        gene_id_to_np =  {}
        db = {}

        for i, row in enumerate(gpff_list):

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
                gene_id = data_value["GeneID"]

                if gene_id not in gene_id_to_np:
                    gene_id_to_np[gene_id] = []

                gene_id_to_np[gene_id].append(np_id)

                data_value.pop("GeneID", None)
            
            db[np_id][data_key] = data_value

        return db, gene_id_to_np

    @staticmethod
    def get_gene_details_from_gpff_list(gene_to_np: (str, list), biomart_helper: BiomartHelper, db):
        gene_id, np_ids = gene_to_np

        row = {}
        row["gene_id"] = [gene_id]
        row["protein_numbers"] = []
        row["gene_synonyms"] = []
        row["gene_symbol"] = []
        row["other_ids"] = []
        row["gene_description"] = []

        for np_id in np_ids:
            protein_data = db[np_id]

            if "gene" in protein_data:
                row["gene_symbol"].append(protein_data["gene"])
                if "locus_tag" in protein_data:
                    row["gene_synonyms"] += protein_data["locus_tag"]
            elif "locus_tag" in protein_data:
                row["gene_symbol"] += protein_data["locus_tag"]

            for key, value in protein_data["db_xrefs"].items():

                if (
                    key in biomart_helper.allowed_db_list
                    and value not in BiomartHelper.allowed_db_list
                ):
                    row["other_ids"].append(value)
                elif key == "EnsemblGenomes-Tr":
                    row["gene_synonyms"].append(value)

            if "gene_synonym" in protein_data:
                row["gene_synonyms"] += protein_data["gene_synonym"]

            row["protein_numbers"] += [np_id]

        biomart_row = biomart_helper.data[
            biomart_helper.data["NCBI gene ID"] == int(gene_id)
        ]
        biomart_row_by_symbol = biomart_helper.data[
            biomart_helper.data["Gene name"] == row["gene_symbol"][0]
        ]

        row["other_ids"] += biomart_row["Gene stable ID"].dropna().tolist()

        row["gene_description"] += (
            biomart_row["Gene description"].dropna().tolist()
            + biomart_row_by_symbol["Gene description"].dropna().tolist()
        )

        if len(row["gene_description"]) > 1:
            row["gene_description"] = row["gene_description"][0:1]

        row = {key: ",".join(list(set(value))) for key, value in row.items()}

        return row

    @staticmethod
    def generate_combined_gpff_dataframe(gpff_dir):
        df_all = pd.DataFrame()
        biomart_helper = BiomartHelper()

        for species_file in os.listdir(gpff_dir):
            if ".list" not in species_file:
                continue

            db, gene_id_to_np = GpffHelper.process_gpff_list(path.join(gpff_dir, species_file))
 
            rows = run_parallel(
                func=GpffHelper.get_gene_details_from_gpff_list,
                args=gene_id_to_np.items(),
                db=db,
                biomart_helper=biomart_helper,
            )

            df = pd.DataFrame(rows)
            df["species_id"] = species_file.replace(".gpff.list", "").replace("_", " ")
            df_all = pd.concat([df_all, df], axis=0)

        return df_all

    @staticmethod
    def save_gpff_to_db(df_all, con):
        cur = con.cursor()

        for index, row in df_all.iterrows():
            gene_id = row["gene_id"]

            for key, values in row.to_dict().items():
                if key == "gene_id":
                    continue
                if key[-1] == "s":
                    key = key[:-1]
                try:
                    if values == "nan":
                        continue

                    if key == "gene_description":
                        cur.execute(
                            "INSERT INTO ncbi_gene_meta(ncbi_gene_id, meta_key, meta_value) VALUES(%s, %s, %s)",
                            (gene_id, key, values),
                        )
                        continue

                    for value in str(values).split(","):
                        if value == "nan":
                            continue
                        cur.execute(
                            "INSERT INTO ncbi_gene_meta(ncbi_gene_id, meta_key, meta_value) VALUES(%s, %s, %s)",
                            (gene_id, key, value),
                        )
                except Exception as e:
                    print(e, gene_id, gene_id, key, value)

        con.commit()

    @staticmethod
    def run_pipeline(staging_path, con):
        biomart_helper = BiomartHelper()

        work_dir = path.join(staging_path, "mapping", "gpff_work")
        list_dir, gpff_dir = (
            path.join(work_dir, "list"),
            path.join(work_dir, "gpff"),
        )

        for p in [list_dir, gpff_dir]:
            pathlib.Path(p).mkdir(parents=True, exist_ok=True)

        GpffHelper.download_files(gpff_dir)
        GpffHelper.gpff_to_list(gpff_dir, list_dir)
        df_gpff = GpffHelper.generate_combined_gpff_dataframe(list_dir)
        
        GpffHelper.save_gpff_to_db(df_gpff, con)

if __name__ == "__main__":

    con = connection_loader()
    config = config_loader()
    GpffHelper.run_pipeline(config["staging_path"], con)
