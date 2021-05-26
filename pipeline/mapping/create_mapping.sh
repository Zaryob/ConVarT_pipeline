#!/bin/bash
BASEDIR=$(dirname "$0")

parse_yaml() {
   local prefix=$2
   local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs=$(echo @|tr @ '\034')
   sed -ne "s|^\($s\)\($w\)$s:$s\"\(.*\)\"$s\$|\1$fs\2$fs\3|p" \
        -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p"  $1 |
   awk -F$fs '{
      indent = length($1)/2;
      vname[indent] = $2;
      for (i in vname) {if (i > indent) {delete vname[i]}}
      if (length($3) > 0) {
         vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
         printf("%s%s%s=\"%s\"\n", "'$prefix'",vn, $2, $3);
      }
   }'
}


creating_protein_mapping_file() {
    #Versin GRCh38.p12:::: zcat GCF_000001405.38_GRCh38.p12_protein.gpff.gz | awk '/REFSEQ: accession/ {print $4}' | awk '$0!="NC_012920.1"' > NM.list
    #Versin GRCh38.p12:::: zcat GCF_000001405.38_GRCh38.p12_protein.gpff.gz | awk '/VERSION/ {print $2}' | awk '$0!~/^YP_/' > NP.list
    zcat $STAGING_PATH"/mapping/gpff_work/gpff/Homo_sapiens.gpff.gz" | awk '/REFSEQ: accession/ {print $4}' | awk -F"." '{print $1}' > $STAGING_PATH"/mapping/NM.list"
    zcat $STAGING_PATH"/mapping/gpff_work/gpff/Homo_sapiens.gpff.gz" | awk '/VERSION/ {print $2}' > $STAGING_PATH"/mapping/NP.list"
    zcat $STAGING_PATH"/mapping/gpff_work/gpff/Homo_sapiens.gpff.gz" | awk '/\/db_xref="GeneID:/ {print $1}' | cut -c18- | tr -d '"' > $STAGING_PATH"/mapping/geneID.list"
    paste $STAGING_PATH"/mapping/NM.list" $STAGING_PATH"/mapping/NP.list" > $STAGING_PATH"/mapping/NM_to_NP.list"
    paste $STAGING_PATH"/mapping/NM_to_NP.list" $STAGING_PATH"/mapping/geneID.list" > $STAGING_PATH"/mapping/NM_NP_GeneID.list"
    rm $STAGING_PATH"/mapping/NP.list" $STAGING_PATH"/mapping/NM.list" $STAGING_PATH"/mapping/NM_to_NP.list" $STAGING_PATH"/mapping/geneID.list"
}

#CREATE THE MAPPING FILE THAT MAPS THE NM_ NUMBERS CORRESPOND WITH NP_ NUMBERS AND GeneIDs
creating_protein_mapping_file() {
    #Versin GRCh38.p12:::: zcat GCF_000001405.38_GRCh38.p12_protein.gpff.gz | awk '/REFSEQ: accession/ {print $4}' | awk '$0!="NC_012920.1"' > NM.list
    #Versin GRCh38.p12:::: zcat GCF_000001405.38_GRCh38.p12_protein.gpff.gz | awk '/VERSION/ {print $2}' | awk '$0!~/^YP_/' > NP.list
    zcat $STAGING_PATH"/mapping/gpff_work/gpff/Homo_sapiens.gpff.gz" | awk '/REFSEQ: accession/ {print $4}' | awk -F"." '{print $1}' > $STAGING_PATH"/mapping/NM.list"
    zcat $STAGING_PATH"/mapping/gpff_work/gpff/Homo_sapiens.gpff.gz" | awk '/VERSION/ {print $2}' > $STAGING_PATH"/mapping/NP.list"
    zcat $STAGING_PATH"/mapping/gpff_work/gpff/Homo_sapiens.gpff.gz" | awk '/\/db_xref="GeneID:/ {print $1}' | cut -c18- | tr -d '"' > $STAGING_PATH"/mapping/geneID.list"
    paste $STAGING_PATH"/mapping/NM.list" $STAGING_PATH"/mapping/NP.list" > $STAGING_PATH"/mapping/NM_to_NP.list"
    paste $STAGING_PATH"/mapping/NM_to_NP.list" $STAGING_PATH"/mapping/geneID.list" > $STAGING_PATH"/mapping/NM_NP_GeneID.list"
    rm $STAGING_PATH"/mapping/NP.list" $STAGING_PATH"/mapping/NM.list" $STAGING_PATH"/mapping/NM_to_NP.list" $STAGING_PATH"/mapping/geneID.list"
}

#NEW CONVERSION TABLE
conversion(){
    declare -A map
    TEXT=$(cat ${PROJECT}/db/mapping/NewCurated_ENSTvsGENEID.csv | cut -d , -f 1)

    while read line;
    do
       map[$line]=1;
    done  <<< "$TEXT"

    IFS=$'\n'; for ENST_number in $(cat "/opt/current_project/db/mapping/WHOLE_ENST_IDS.csv"); 
    do
        #count_in_exist_file=`cat $PROJECT"/db/mapping/NewCurated_ENSTvsGENEID.csv" | awk -v FS="," '$1=="'$ENST_number'"' |wc -l`
        #count_in_exist_file=`echo $TEXT | grep "${ENST_number},"`
        if [[ -v "map[$ENST_number]" ]] ; then
          echo "$ENST_number already exists";
          continue;
        else
            echo "$ENST_number is being added";
        fi
        
        hgnc_id=`cat $STAGING_PATH"/mapping/ENST_ENSTV_GeneID.list" | awk -F"," '$1=="'$ENST_number'" {print $4}' | sort -u`
        hgnc_id_number=`cat $STAGING_PATH"/mapping/ENST_ENSTV_GeneID.list" | awk -F"," '$1=="'$ENST_number'" {print $4}' | sort -u | wc -l`
        if [ ! -z "$hgnc_id" ] && [ "$hgnc_id_number" == "1" ]; then
            term="HGNC:$hgnc_id"
            gene_id=`cat $STAGING_PATH"/mapping/HGNC.list" | awk -F"\t" '$1=="'$term'" {print $5}' | sort -u | head -n 1`
        else
            missing_gene_name=`cat $STAGING_PATH"/mapping/ENST_ENSTV_GeneID.list" | awk -F"," '$1=="'$ENST_number'" {print $3}' | sort -u`
            if [ -z "$missing_gene_name" ]; then
                continue;
            fi
            
            esearch -db gene -query "($missing_gene_name[Gene Name]) AND 9606[Taxonomy ID]" | efetch -format tabular > $TMP_PATH"/entrez.csv"
            gene_id=`tail -n +2 $TMP_PATH"/entrez.csv" | awk -F"\t" '$6=="'$missing_gene_name'" {print $3}' | sort -u`
            if [ -z "$gene_id" ]; then
                esearch -db gene -query "($ENST_number) AND 9606[Taxonomy ID]" | efetch -format tabular > $TMP_PATH"/entrez.csv"
                gene_id=`tail -n +2 $TMP_PATH"/entrez.csv" | awk -F"\t" '{print $3}' | sort -u`
            fi

        fi
        echo $ENST_number,$gene_id >> $STAGING_PATH"/mapping/NewCurated_ENSTvsGENEID.csv"
    done
    rm  $TMP_PATH"/entrez.csv"
}

creating_protein_mapping_file()
conversion()
