#!/bin/bash

#PROJECT ENVIRONMENT
PROJECT="/opt/current_project"
DATABASE="/opt/current_project/db"

#CLINVAR
retrieve_and_process_ClinVar() {

    #WHOLE DATA SET {329733 rows}
    #├── GENEID == -1 {751 rows} (no problem with NM_ which means all the rows have NM numbers)
    #    ├── "NT expansion" exist {3 rows}  /do it manuel/
    #    ├── NT expansion NOT exist {748 rows}
    #├── GENEID != -1 {328982 rows}
    #├──NT expansion && NM_  exist {21 rows} /do it manuel/
    #├──NT expansion && NM_ NOT Exist {52 rows} /do it manuel/
    #├──NT expansion yok && NM_ exist {328909 rows}


    #downloaded_time=`date | awk -F" " '{print $3""$2""$6}'`
    #echo -e "________________________________________________________________________________________ \n New ClinVar file is downloading right now... \n____________________________________________________________________________________ \n \n"
    #wget "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
    #mv "variant_summary.txt.gz" $DATABASE"/clinvar/variant_summary_"$downloaded_time".txt.gz"
    touch $DATABASE"/clinvar/firstFormat4Şub2019.csv" 
    zcat $DATABASE"/clinvar/variant_summary_4Şub2019.txt.gz" | awk -v FS="\t" -v OFS="\t"  '{print $4,$1,$5,$10,$12,$31,$2,$3,$7,$9,$14,$24,$25}' | head -n 1 > $DATABASE"/clinvar/firstFormat4Şub2019.csv"
    zcat $DATABASE"/clinvar/variant_summary_4Şub2019.txt.gz" | awk -v FS="\t" -v OFS="\t"  '$1!="-1" && $17=="GRCh37" {print $4,$1,$5,"rs"$10,$12,$31,$2,$3,$7,$9,$14,$24,$25}' | awk '/\(p./' >> $DATABASE"/clinvar/firstFormat4Şub2019.csv"

    #Create Folders
    mkdir $DATABASE"/clinvar/gene_id_no/"  
    mkdir $DATABASE"/clinvar/gene_id_yes/"  

    #GeneID_NOTExist && NT expansion_NOTExist {748}
    touch $DATABASE"/clinvar/gene_id_no/geneIDno_NTno.txt"
    tail -n +2 $DATABASE"/clinvar/firstFormat4Şub2019.csv" | awk -v FS="\t" -v OFS="\t" '$1=="-1" && $7!="NT expansion"' > $DATABASE"/clinvar/gene_id_no/geneIDno_NTno.txt"

    #GeneID_NOTExist && NT expansion_NOTExist {3}
    touch $DATABASE"/clinvar/gene_id_no/geneIDno_NTyes.txt"
    tail -n +2 $DATABASE"/clinvar/firstFormat4Şub2019.csv" | awk -v FS="\t" -v OFS="\t" '$1=="-1" && $7=="NT expansion"' > $DATABASE"/clinvar/gene_id_no/geneIDno_NTyes.txt"

    #GeneID_Exist && NT expansion_Exist && NM_Exist {21}
    touch $DATABASE"/clinvar/gene_id_yes/geneIDyes_NTyes_NMyes.txt"
    tail -n +2 $DATABASE"/clinvar/firstFormat4Şub2019.csv" | awk -v FS="\t" -v OFS="\t" '$1!="-1" && $7=="NT expansion"' > $DATABASE"/clinvar/gene_id_yes/geneIDyes_NTyes_NMyes.txt"

    #GeneID_Exist && NT expansion_NOTExist && NM_NOTExist {52}
    touch $DATABASE"/clinvar/gene_id_yes/fullNamesForNTno_tmp.txt" 
    touch $DATABASE"/clinvar/gene_id_yes/NMyes_tmp.txt" 
    tail -n +2 $DATABASE"/clinvar/firstFormat4Şub2019.csv" | awk -v FS="\t" -v OFS="\t" '$1!="-1" && $7!="NT expansion"' > $DATABASE"/clinvar/gene_id_yes/fullNamesForNTno_tmp.txt"
    tail -n +2 $DATABASE"/clinvar/firstFormat4Şub2019.csv" | awk '/NM_/' | awk -v FS="\t" -v OFS="\t" '$1!="-1" && $7!="NT expansion"' > $DATABASE"/clinvar/gene_id_yes/NMyes_tmp.txt"
    awk 'FNR==NR{a[$0]; next} !($0 in a)' $DATABASE"/clinvar/gene_id_yes/NMyes_tmp.txt" $DATABASE"/clinvar/gene_id_yes/fullNamesForNTno_tmp.txt" > $DATABASE"/clinvar/gene_id_yes/geneIDyes_NTno_NMno.txt"
    rm $DATABASE"/clinvar/gene_id_yes/NMyes_tmp.txt" $DATABASE"/clinvar/gene_id_yes/fullNamesForNTno_tmp.txt"

    #GeneID_Exist && NT expansion_NOTExist && NM_Exist {328909}
    touch $DATABASE"/clinvar/gene_id_yes/geneIDyes_NTno_NMyes.txt"
    tail -n +2 $DATABASE"/clinvar/firstFormat4Şub2019.csv" | awk '/NM_/' | awk -v FS="\t" -v OFS="\t" '$1!="-1" && $7!="NT expansion"' > $DATABASE"/clinvar/gene_id_yes/geneIDyes_NTno_NMyes.txt"

    #PARSING PROTEIN CHANGES
    touch $DATABASE"/clinvar/gene_id_yes/noSpaceBeforeProteinAlt.txt"
    touch $DATABASE"/clinvar/gene_id_yes/VariationIDsOfSpaced.txt"
    tail -n +2 $DATABASE"/clinvar/firstFormat4Şub2019.csv" | awk -v FS="\t" '{print $8}' | awk -v FS=" " '$2=="" {print $0}' > $DATABASE"/clinvar/gene_id_yes/noSpaceBeforeProteinAlt.txt"
    cat $DATABASE"/clinvar/gene_id_yes/noSpaceBeforeProteinAlt.txt" | while read name
    do
        awk -v FS="\t" '$8=="'$name'" {print $6}' >> $DATABASE"/clinvar/gene_id_yes/VariationIDsOfSpaced.txt"
    done

    #PARSING PROTEIN CHANGES
    tail -n +2 $DATABASE"/clinvar/gene_id_yes/geneIDyes_NTno_NMyes.txt" | awk -v FS="\t" '{print $8}' | awk -v FS=" " '$2=="" {print $0}' > $DATABASE"/clinvar/gene_id_yes/noSpaceBeforeProteinAlt.txt"
    cat $DATABASE"/clinvar/gene_id_yes/noSpaceBeforeProteinAlt.txt" | while read name
    do
        cat $DATABASE"/clinvar/firstFormat4Şub2019.csv" | awk -v FS="\t" '$8=="'$name'" {print $6}' >> $DATABASE"/clinvar/gene_id_yes/VariationIDsOfSpaced.txt"
    done

    cat $DATABASE"/clinvar/gene_id_yes/VariationIDsOfSpaced.txt" | while read variationID
    do
        # in order to retrieve the data which has not space before protein alteration on geneIDyes_NTno_NMyes.txt {7 rows}
        cat $DATABASE"/clinvar/firstFormat4Şub2019.csv" | awk -v FS="\t" '$6=="'$variationID'" {print $0}' >> $DATABASE"/clinvar/gene_id_yes/spaceProblemForParserFullRecords.txt" 
    done

    # to be able to retrieve and overwrite geneIDyes_NTno_NMyes.txt with spaces (before protein alteration) {328909 -> 328909-7 ==> 328902 rows}
    cat $DATABASE"/clinvar/gene_id_yes/geneIDyes_NTno_NMyes.txt" | awk -v FS="\t" '$6!="21087" && $6!="24813" && $6!="217159" && $6!="218831" && $6!="424933" && $6!="440232" && $6!="498619" {print $0}' > $DATABASE"/clinvar/gene_id_yes/temporaryForNew_geneIDyes_NTno_NMyes.txt" 

    rm $DATABASE"/clinvar/gene_id_yes/noSpaceBeforeProteinAlt.txt" 
    rm $DATABASE"/clinvar/gene_id_yes/VariationIDsOfSpaced.txt"
    rm $DATABASE"/clinvar/gene_id_yes/geneIDyes_NTno_NMyes.txt"
    mv $DATABASE"/clinvar/gene_id_yes/temporaryForNew_geneIDyes_NTno_NMyes.txt" $DATABASE"/clinvar/gene_id_yes/geneIDyes_NTno_NMyes.txt" 


    #Split the numbers
    cat $DATABASE"/clinvar/gene_id_yes/geneIDyes_NTno_NMyes.txt" | awk -F"\t" '{print $8}' | awk -F" " '{print $2}' >  $DATABASE"/clinvar/gene_id_yes_variations.txt"
    cat $DATABASE"/clinvar/gene_id_yes/geneIDyes_NTno_NMyes.txt" | awk -F"\t" '{print $8}' | awk -F" " '{print $2}' | sed 's/[0-9][0-9]*/ & /' | awk -v FS=" " '{print $2}' > $DATABASE"/clinvar/gene_id_yes_numbers.txt"
    paste $DATABASE"/clinvar/gene_id_yes/geneIDyes_NTno_NMyes.txt" $DATABASE"/clinvar/gene_id_yes_numbers.txt" > $DATABASE"/clinvar/gene_id_yes/withPositions.txt"
    paste $DATABASE"/clinvar/gene_id_yes/withPositions.txt" $DATABASE"/clinvar/gene_id_yes_variations.txt" > $DATABASE"/clinvar/final/geneIDyes_NTno_NMyes.csv"
    rm $DATABASE"/clinvar/gene_id_yes/withPositions.txt" $DATABASE"/clinvar/gene_id_yes_numbers.txt" $DATABASE"/clinvar/gene_id_yes_variations.txt"
}

#CREATE THE MAPPING FILE THAT MAPS THE NM_ NUMBERS CORRESPOND WITH NP_ NUMBERS AND GeneIDs
creating_protein_mapping_file() {
    #Versin GRCh38.p12:::: zcat GCF_000001405.38_GRCh38.p12_protein.gpff.gz | awk '/REFSEQ: accession/ {print $4}' | awk '$0!="NC_012920.1"' > NM.list
    #Versin GRCh38.p12:::: zcat GCF_000001405.38_GRCh38.p12_protein.gpff.gz | awk '/VERSION/ {print $2}' | awk '$0!~/^YP_/' > NP.list
    zcat $DATABASE"/gpff/Homo_sapiens.gpff.gz" | awk '/REFSEQ: accession/ {print $4}' | awk -F"." '{print $1}' > $DATABASE"/mapping/NM.list"
    zcat $DATABASE"/gpff/Homo_sapiens.gpff.gz" | awk '/VERSION/ {print $2}' > $DATABASE"/mapping/NP.list"
    zcat $DATABASE"/gpff/Homo_sapiens.gpff.gz" | awk '/\/db_xref="GeneID:/ {print $1}' | cut -c18- | tr -d '"' > $DATABASE"/mapping/geneID.list"
    paste $DATABASE"/mapping/NM.list" $DATABASE"/mapping/NP.list" > $DATABASE"/mapping/NM_to_NP.list"
    paste $DATABASE"/mapping/NM_to_NP.list" $DATABASE"/mapping/geneID.list" > $DATABASE"/mapping/NM_NP_GeneID.list"
    rm $DATABASE"/mapping/NP.list" $DATABASE"/mapping/NM.list" $DATABASE"/mapping/NM_to_NP.list" $DATABASE"/mapping/geneID.list"
}

#COMPARE SEQUENCE FROM ClinVar and gnomAD
compare_seqs() {
    echo "seq_from_clinvar,geneID,seq_from_gnomAD,needle_score,stretcher_score" > $DATABASE"/mapping/SequenceMatch_ClinVar_gnomAD.list"
    ls -l /opt/current_project/results/clinvar_seqs/*_ClinVar.fasta | awk -F"/" '{print $6}' | awk -F"_ClinVar.fasta" '{print $1}' | while read NP_number
    do
        gene_id=`cat $DATABASE"/mapping/NM_NP_GeneID.list" | awk -v FS="\t" '$2=="'$NP_number'" {print $3}'`
        if [ -z "$gene_id" ]; then
            echo $NP_number",noGeneID,-,0,0" >> $DATABASE"/mapping/SequenceMatch_ClinVar_gnomAD.list"
        else
            cat $DATABASE"/mapping/NewCurated_ENSTvsGENEID.csv" | awk -v FS="," '$2=="'$gene_id'" {print $1}' | while read ENST_number
            do
                if [ -z "$ENST_number" ]; then
                    echo $NP_number","$gene_id",noENSTnumber,0,0"  >> $DATABASE"/mapping/SequenceMatch_ClinVar_gnomAD.list"
                else
                    if [ -e "$PROJECT/results/gnomad_seqs/"$ENST_number"_gnomAD.fasta" ]; then
                    #If gnomad_seqs folder consists of this enst_number, continue. Otherwise, pass.
                    needle -asequence $PROJECT"/results/clinvar_seqs/"$NP_number"_ClinVar.fasta" -bsequence $PROJECT"/results/gnomad_seqs/"$ENST_number"_gnomAD.fasta" -gapopen 10.0 -gapextend 0.5 -endopen 10.0 -endextend 0.5 -sprotein1 -sprotein2 -aformat3 pair -auto -stdout -outfile $PROJECT"/tmp/needle.txt"
                    needle_score=`cat $PROJECT"/tmp/needle.txt" | awk -F[:\(\)] '/# Identity:/  {print $3}'`
                    stretcher -asequence $PROJECT"/results/clinvar_seqs/"$NP_number"_ClinVar.fasta" -bsequence $PROJECT"/results/gnomad_seqs/"$ENST_number"_gnomAD.fasta" -gapopen 12 -gapextend 2 -sprotein1 -sprotein2 -aformat3 pair -auto -stdout -outfile $PROJECT"/tmp/stretcher.txt"
                    stretcher_score=`cat $PROJECT"/tmp/stretcher.txt" | awk -F[:\(\)] '/# Identity:/  {print $3}'`
                    echo $NP_number","$gene_id","$ENST_number","$needle_score","$stretcher_score >> $DATABASE"/mapping/SequenceMatch_ClinVar_gnomAD.list"
                    fi
                fi
            done
        fi
    done
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
        
        hgnc_id=`cat $DATABASE"/mapping/ENST_ENSTV_GeneID.list" | awk -F"," '$1=="'$ENST_number'" {print $4}' | sort -u`
        hgnc_id_number=`cat $DATABASE"/mapping/ENST_ENSTV_GeneID.list" | awk -F"," '$1=="'$ENST_number'" {print $4}' | sort -u | wc -l`
        if [ ! -z "$hgnc_id" ] && [ "$hgnc_id_number" == "1" ]; then
            term="HGNC:$hgnc_id"
            gene_id=`cat $DATABASE"/mapping/HGNC.list" | awk -F"\t" '$1=="'$term'" {print $5}' | sort -u | head -n 1`
        else
            missing_gene_name=`cat $DATABASE"/mapping/ENST_ENSTV_GeneID.list" | awk -F"," '$1=="'$ENST_number'" {print $3}' | sort -u`
            if [ -z "$missing_gene_name" ]; then
                continue;
            fi
            
            esearch -db gene -query "($missing_gene_name[Gene Name]) AND 9606[Taxonomy ID]" | efetch -format tabular > $PROJECT"/tmp/entrez.csv"
            gene_id=`tail -n +2 $PROJECT"/tmp/entrez.csv" | awk -F"\t" '$6=="'$missing_gene_name'" {print $3}' | sort -u`
            if [ -z "$gene_id" ]; then
                esearch -db gene -query "($ENST_number) AND 9606[Taxonomy ID]" | efetch -format tabular > $PROJECT"/tmp/entrez.csv"
                gene_id=`tail -n +2 $PROJECT"/tmp/entrez.csv" | awk -F"\t" '{print $3}' | sort -u`
            fi

        fi
        echo $ENST_number,$gene_id >> $PROJECT"/db/mapping/NewCurated_ENSTvsGENEID.csv"
    done
    rm  $PROJECT"/tmp/entrez.csv"
}

#Create CONVART Human Transcript File and BLAST_DB
create_proteome() {
    if [ -f "$DATABASE/blast_db/convart_curated_fasta.fasta.fasta" ]; then rm -r $DATABASE"/blast_db/"; mkdir $DATABASE"/blast_db/"; fi
    #cat $DATABASE"/domains/domain_gnomaAD.fasta" $DATABASE"/domains/domain_ClinVar.fasta" >> $DATABASE"/blast_db/convartProteome.fasta"
    makeblastdb -in "$DATABASE/blast_db/convart_curated_fasta.fasta" -out "$DATABASE/blast_db/convartProteome" -title "$DATABASE/blast_db/convartProteome" -dbtype prot 
}

#Generate Domains
generate_domains() {
    /opt/CilioGenics/domain/PfamScan/pfam_scan.pl -dir /opt/CilioGenics/domain/db/ -fasta $DATABASE"/blast_db/convart_curated_fasta.fasta" > $DATABASE"/domains/Raw_Domains.list"
    tail -n +29 $DATABASE"/domains/Raw_Domains.list" | awk -v FS=" " -v OFS="\t" '{print $1,$6,$7,$8,$2,$3,$12,$13,$15}' > $DATABASE"/domains/Domains.list"
}
