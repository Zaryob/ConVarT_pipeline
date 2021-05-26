import sys
import os
from envyaml import EnvYAML
import pymysql
import pandas as pd
from pybiomart import Server
import subprocess

def config_loader():
    dir_path = os.path.dirname(os.path.realpath(__file__))
    config_file = os.path.join(dir_path, "..", "config", "commons.yml")
    config = EnvYAML(config_file)

    return config

def connection_loader():

    config = config_loader()
    db_config = config["database"]

    con = pymysql.connect(
        host=db_config["host"],
        user=db_config["username"],
        passwd=db_config["password"],
        db=db_config["database"]
    )

    return con

def catched_subprocess(command):
    try:
        return subprocess.check_output(command, shell=True)
    except Exception as e:
        print("Error", e, command)
        return ''.encode()

class BiomartHelper:
    species_to_dataset = pd.DataFrame(
        [
            {
                "species_name": "Mus_musculus",
                "table": "mapping_mouse",
                "dataset": "mmusculus_gene_ensembl",
            },
            {
                "species_name": "Rattus_norvegicus",
                "table": "mapping_rat",
                "dataset": "rnorvegicus_gene_ensembl",
            },
            {
                "species_name": "Danio_rerio",
                "table": "mapping_zebrafish",
                "dataset": "drerio_gene_ensembl",
            },
            {
                "species_name": "Homo_sapiens",
                "table": "mapping_human",
                "dataset": "hsapiens_gene_ensembl",
            },
            {
                "species_name": "Pan_troglodytes",
                "table": "mapping_chimp",
                "dataset": "ptroglodytes_gene_ensembl",
            },
            {
                "species_name": "Drosophila_melanogaster",
                "table": "mapping_fruitfly",
                "dataset": "dmelanogaster_gene_ensembl",
            },
            {
                "species_name": "Macaca_mulatta",
                "table": "mapping_macaque",
                "dataset": "mmulatta_gene_ensembl",
            },
            {
                "species_name": "Xenopus_tropicalis",
                "table": "mapping_frog",
                "dataset": "xtropicalis_gene_ensembl",
            },
            {
                "species_name": "Caenorhabditis_elegans",
                "table": "mapping_worm",
                "dataset": "celegans_gene_ensembl",
            },
        ]
    )
    allowed_db_list = [
        "Xenbase",
        "Wormbase",
        "HGNC",
        "FLYBASE",
        "EnsembleGenomes-Gn",
        "MGI",
        "ZFIN",
    ]

    data = None

    def __init__(self, species_to_dataset=None):
        if species_to_dataset is not None:
            self.species_to_dataset = species_to_dataset

        self.fetch_data()

    def fetch_data(self):
        data = pd.DataFrame()

        server = Server(host="http://www.ensembl.org")

        datasets = self.species_to_dataset["dataset"].tolist()
        for dataset_name in datasets:
            dataset = server.marts["ENSEMBL_MART_ENSEMBL"].datasets[dataset_name]
            dataset_df = dataset.query(
                attributes=[
                    "ensembl_gene_id",
                    "entrezgene_id",
                    "description",
                    "external_gene_name",
                ]
            )

            data = pd.concat([data, dataset_df], axis=0)

        self.data = data

