#!/bin/bash

BASEDIR=$(dirname "$0")
source ~/.bash_profile
source ~/.bashrc
echo $VEP_PATH;
eval $(sed -E 's/:( |\t){1,}(.*?)$/=\2;\n/g' $BASEDIR/../config.yml >&1)

perl5.26.1 ~/vep/vep -i $PROJECT/dbSNP/common_all_20180418.vcf  -o $PROJECT/dbSNP/result.txt \
                --cache -offline -synonym ~/.vep/homo_sapiens/97_GRCh37/chr_synonyms.txt --hgvs --fasta \
                $VEP_DATA/homo_sapiens/$VER\_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
                --no_check_variants_order
