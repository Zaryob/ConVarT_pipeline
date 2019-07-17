#!/bin/bash

<< ///
	   _____                          _     _____           _           _
	  / ____|                        | |   |  __ \         (_)         | |
	 | |    _   _ _ __ _ __ ___ _ __ | |_  | |__) | __ ___  _  ___  ___| |_
	 | |   | | | | '__| '__/ _ \ '_ \| __| |  ___/ '__/ _ \| |/ _ \/ __| __|
	 | |___| |_| | |  | | |  __/ | | | |_  | |   | | | (_) | |  __/ (__| |_
	  \_____\__,_|_|  |_|  \___|_| |_|\__| |_|   |_|  \___/| |\___|\___|\__|
        	                                              _/ |
                	                                     |__/
	Current Project v0.1
	____________________

	Created Date: 01.02.2019
	Last Update: 01.02.2019
		* First commit
///

#PROJECT ENVIRONMENT
PROJECT="/opt/current_project"
DATABASE="/opt/current_project/db"



retrieve_and_process_ClinVar() {


	#WHOLE DATA SET {329733 rows}
	#├── GENEID == -1 {751 rows} (no problem with NM_ which means all the rows have NM numbers)
	#    ├── "NT expansion" exist {3 rows}  /do it manuel/
	#	 ├── NT expansion NOT exist {748 rows}
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
	zcat $DATABASE"/clinvar/variant_summary_4Şub2019.txt.gz" | awk -v FS="\t" -v OFS="\t"  '$1!="-1" && $17=="GRCh38" {print $4,$1,$5,"rs"$10,$12,$31,$2,$3,$7,$9,$14,$24,$25}' | awk '/\(p./' >> $DATABASE"/clinvar/firstFormat4Şub2019.csv"

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
    
	#preProcesses for Parsing ClinVar    
	#GECE IPTALI!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!: touch $DATABASE"/clinvar/gene_id_yes/noSpaceBeforeProteinAlt.txt"
	#GECE IPTALI!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!: touch $DATABASE"/clinvar/gene_id_yes/VariationIDsOfSpaced.txt"

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

retrieve_and_process_ClinVar

#DETECT THE ORGANISMS IN THE FOLDER
detect_organisms() {
        ls -l $DATABASE"/pep-fasta/" | awk '{print $9}' | tail -n +2 > $DATABASE"/organisms.txt"
}

#CREATE BLAST INDEX
create_blast_index() {
	cat $DATABASE"/organisms.txt" | while read organism
	do
		organism_name=$(echo $organism | awk -F"." '{print $1}')
		makeblastdb -in $DATABASE"/pep-fasta/"$organism  -out $DATABASE"/pep-blast/"$organism_name -title $DATABASE"/pep-blast/"$organism_name -dbtype prot
	done
}

<< CODE_DOWN
#PREPARE QUERY FILES FROM FASTA FILE BY SELECTING THE LONGEST GENES
prepare_query_files() {
	cat $DATABASE"/organisms.txt" | while read organism
	do
		organism_name=$(echo $organism | awk -F"." '{print $1}')
		awk -f $PROJECT"/pipeline/reducegene.awk" $DATABASE"/pep-fasta/"$organism  > $DATABASE"/pep-query/"$organism_name
	done
}
CODE_DOWN

#RETRIEVE THE SEQUENCES OF THE PROTEINS ACROSS MODEL ORGANISMS
sequence_retrieve() {
	cat $DATABASE"/others/selectedgenes.txt" | while read gene_name
	do
		gene_symbol=$(echo $gene_name)
		searchterm="gene_symbol:$gene_symbol"
		awk -v RS=">" '$8=="'$searchterm'" {print ">"$0}' $DATABASE"/pep-query/human.fa" > $PROJECT"/results/longest_products_from_human.fa"
	done
}

#TIME TO BLASTS
blast_hit() {
	cat $DATABASE"/others/organisms.txt" | while read organism
	organism_name=$(echo $organism | awk -F"." '{print $1}')
	do
	echo -e "_________________________________________________ \n$organism_name >> $organism_name _____________________________________________________________"
	cat $PROJECT"/results/longest_products_from_human.fa" | parallel --blocksize 10K --memfree 2G --recstart '>' --pipe \
	blastp -num_threads 16 -outfmt 6 -db $DATABASE"/pep-blast/"$organism_name  -query - -max_target_seqs 5 > $PROJECT"/results/blast_hits/"$organism_name".csv"
	done
}