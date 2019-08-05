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

///

#PROJECT ENVIRONMENT
PROJECT="/opt/current_project"
DATABASE="/opt/current_project/db"

#CLINVAR
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

#RETRIEVE CLINVAR SEQUENCE
ret_clinvar_seq_domains() {
	cat $DATABASE"/clinvar/final/ClinVar.csv" | awk -F"\t" '{print $8}' | awk -F"(" '{print $1}' | awk '$0~/^NM_/' | awk -F"." '{print $1}' | sort -u > $DATABASE"/mapping/ClinVar_NM_unique.list"

	cat $DATABASE"/mapping/ClinVar_NM_unique.list" | while read ClinVar_NM_number
	do
		line=`cat $DATABASE"/mapping/NM_NP_GeneID.list" | awk -F"\t" '$1=="'$ClinVar_NM_number'" {print $0}'`
		NP_number=`echo $line | awk '{print $2}'`
		gene_id=`echo $line | awk '{print $3}'`
		seq=`zcat $DATABASE"/proteins/Homo_sapiens.faa.gz" | awk -v RS="(^|\n)>" '$1=="'$NP_number'" {print ">"$0}'`
		if [ -z "$NP_number" ]; then
			echo $NP_number "is NOT FOUND for " $ClinVar_NM_number >> "logs/ret_clinvar_seq.log"
		else
			echo $ClinVar_NM_number "--> GeneID: " $gene_id "-->" $NP_number "[for Human ClinVar]" >> "logs/ret_clinvar_seq.log"
		fi
		if [ -z "$seq" ]; then
			wget -q "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=protein&report=fasta&id="$NP_number"&extrafeat=null&conwithfeat=on" -O $PROJECT"/results/clinvar_seqs/"$NP_number"_ClinVar.fasta"
			echo "Due to missing in tables or ref. sequence, the file was retrieved from NCBI. <--" $NP_number " <-- GeneID: " $gene_id " <-- " $ClinVar_NM_number >> "logs/ret_clinvar_seq.log"
			sleep 3
		else
			zcat $DATABASE"/proteins/Homo_sapiens.faa.gz" | awk -v RS="(^|\n)>" '$1=="'$NP_number'" {print ">"$0}' > $PROJECT"/results/clinvar_seqs/"$NP_number"_ClinVar.fasta"
		fi
	done

	#DOMAIN FIND
	if [ -f "$DATABASE/domains/domain_ClinVar.fasta"  ]; then rm $DATABASE"/domains/domain_ClinVar.fasta"; fi
	cat $PROJECT"/results/clinvar_seqs/"*"_ClinVar.fasta" > $DATABASE"/domains/domain_ClinVar.fasta"
	if [ -f "$DATABASE/domains/Raw/ClinVar_domains.list"  ]; then rm $DATABASE"/domains/Raw/ClinVar_domains.list"; fi
	/opt/CilioGenics/domain/PfamScan/pfam_scan.pl -dir /opt/CilioGenics/domain/db/ -fasta $DATABASE"/domains/domain_ClinVar.fasta" > $DATABASE"/domains/Raw/ClinVar_domains.list"
	tail -n +29 $DATABASE"/domains/Raw/ClinVar_domains.list" | awk -v FS=" " -v OFS="\t" '{print $1,$6,$7,$8,$2,$3,$12,$13,$15}' > $DATABASE"/domains/ClinVar_Domains.list"
}

#RETRIEVE gnomAD SEQUENCES
ret_gnomad_seq_domains() {
	#DELETE OLDEST FILE
    if [ -f "$PROJECT/pipeline/logs/ret_gnomad_seq.log" ]; then rm $PROJECT"/pipeline/logs/ret_gnomad_seq.log"; fi

    #ALTERNATIVE: zcat Homo_sapiens_Ensembl.fa.gz| awk -v RS="(^|\n)>" '$5=="transcript:ENST00000483390.2" {print ">"$0}' | awk -F" " '/^>/ {$0=$0 " [Homo sapiens]"}1'
    cat $DATABASE"/mapping/gnomAD_ENST_unique.list" | while read ENST_number
    do
            ENST_number_with_version=`cat $DATABASE"/mapping/ENST_ENSTV_GeneID.list" | awk -v FS="," -v OFS="," '$1=="'$ENST_number'" {print $2}' | sort -nk1,1 | head -n 1`
            gene_id=`cat $DATABASE"/mapping/NewCurated_ENSTvsGENEID.csv" | awk -v FS="," -v OFS="," '$1=="'$ENST_number'" {print $2}' | sort -nk1,1 | head -n 1`
            term="transcript:$ENST_number_with_version"
            seq=`zcat $DATABASE"/proteins/Homo_sapiens_Ensembl.fa.gz" | awk -v RS="(^|\n)>" '$5=="'$term'" {print ">"$0}'`
            if [ -z "$ENST_number" ]; then
            	echo "ENST Number (with version) is NOT FOUND for " $ENST_number >> "logs/ret_gnomad_seq.log"
            else
            	echo $ENST_number "--> GeneID: " $gene_id "-->" $ENST_number_with_version "[for Human gnomAD]" >> "logs/ret_gnomad_seq.log"
            fi
            if [ -z "$seq" ]; then
            	wget -q --header='Content-type:text/x-fasta' "https://rest.ensembl.org/sequence/id/"$ENST_number"?type=protein" -O $PROJECT"/results/gnomad_seqs/"$ENST_number"_gnomAD_temp.fasta"
            	cat $PROJECT"/results/gnomad_seqs/"$ENST_number"_gnomAD_temp.fasta" | sed "s/>.*/>$ENST_number_with_version gnomAD:$ENST_number GeneID:$gene_id [Homo sapiens]/" > $PROJECT"/results/gnomad_seqs/"$ENST_number"_gnomAD.fasta"
            	rm $PROJECT"/results/gnomad_seqs/"$ENST_number"_gnomAD_temp.fasta"
				echo "Due to missing in tables or ref. sequence, the file was retrieved from ENSEMBL-Rest-Service. <--" $ENST_number " <-- GeneID: " $gene_id  >> "logs/ret_gnomad_seq.log"
				sleep 3
            else
            	zcat $DATABASE"/proteins/Homo_sapiens_Ensembl.fa.gz" | awk -v RS="(^|\n)>" '$5=="'$term'" {print ">"$0}' | sed "s/>.*/>$ENST_number_with_version gnomAD:$ENST_number GeneID:$gene_id [Homo sapiens]/" > $PROJECT"/results/gnomad_seqs/"$ENST_number"_gnomAD.fasta"
            fi
    done

    #DOMAIN FIND
	if [ -f "$DATABASE/domains/domain_gnomaAD.fasta" ]; then rm $DATABASE"/domains/domain_gnomaAD.fasta"; fi
	#FAILED DUE TO THE HUGE NUMBER OF FILE =>	cat $PROJECT"/results/gnomad_seqs/"*"_gnomAD.fasta" > $DATABASE"/domains/domain_gnomaAD.fasta"
	#FAILED OWING TO AN UNKNOWN REASON =>  		IFS="\n"; for file in find -type f -name "$PROJECT/results/gnomad_seqs/"*"_gnomAD.fasta"; do cat "$file" >> "$DATABASE/domains/domain_gnomaAD.fasta"; done
	list_of_gnomad_files=($(ls $PROJECT"/results/gnomad_seqs/"))
	for gnomad_file in ${list_of_gnomad_files[@]}; do cat $PROJECT"/results/gnomad_seqs/"$gnomad_file >> "$DATABASE/domains/domain_gnomaAD.fasta" ; done
	if [ -f "$DATABASE/domains/Raw/gnomAD_domains.list" ]; then rm $DATABASE"/domains/Raw/gnomAD_domains.list"; fi
	/opt/CilioGenics/domain/PfamScan/pfam_scan.pl -dir /opt/CilioGenics/domain/db/ -fasta $DATABASE"/domains/domain_gnomaAD.fasta" > $DATABASE"/domains/Raw/gnomAD_domains.list"
	tail -n +29 $DATABASE"/domains/Raw/gnomAD_domains.list" | awk -v FS=" " -v OFS="\t" '{print $1,$6,$7,$8,$2,$3,$12,$13,$15}' > $DATABASE"/domains/gnomAD_Domains.list"
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

#GET HOMOLOGY SEQ FOR ClinVar
get_homology_seq_clinvar() {
	if [ -d "$PROJECT/results/seqs" ]; then rm -r $PROJECT"/results/seqs"; fi
	cp $PROJECT"/results/clinvar_seqs" $PROJECT"/results/seqs" -r 
	if [ -f "logs/get_homology_seq_clinvar.log" ]; then rm "logs/get_homology_seq_clinvar.log"; fi

	ls -l /opt/current_project/results/clinvar_seqs/*_ClinVar.fasta | awk -F"/" '{print $6}' | awk -F"_ClinVar.fasta" '{print $1}' | while read NP_number
	do
		gene_id=`cat $DATABASE"/mapping/NM_NP_GeneID.list" | awk -F"\t" '$2=="'$NP_number'" {print $3}'`
		if [ -z "$NP_number" ] || [ -z "$gene_id" ]; then
			echo "NP or GeneID is empty!";
		else
			homology_line=`cat $DATABASE"/orthology/Homology.list" | awk -v FS="\t" -v OFS=";" '$1=="'$gene_id'" {print $1,$2,$3,$4,$5,$6,$7,$8,$9}'`
			#Chimp
			echo $homology_line | awk -v FS=";" '{print $2}' | tr ',' '\n' | while read id
				do
					prot=`cat $DATABASE"/mapping/NP_tables/Pan_troglodytes.list" | awk -v FS="\t" '$6=="'$id'" {print $6,$8,$9}' | sort -rnk3,3 | awk 'NR==1 {print $2}'`
					if [ ! -z "$prot" ]; then
						zcat $DATABASE"/proteins/Pan_troglodytes.faa.gz" | awk -v RS="(^|\n)>" '$1=="'$prot'" {print ">"$0}' >>  $PROJECT"/results/seqs/"$NP_number"_ClinVar.fasta"
						echo $NP_number "-->" $gene_id "-->" $id "-->" $prot "[for Chimp]" >> "logs/get_homology_seq_clinvar.log"
					fi
				done
			#Macaque
			echo $homology_line | awk -v FS=";" '{print $3}' | tr ',' '\n' | while read id
	        do
	            prot=`cat $DATABASE"/mapping/NP_tables/Macaca_mulatta.list" | awk -v FS="\t" '$6=="'$id'" {print $6,$8,$9}' | sort -rnk3,3 | awk 'NR==1 {print $2}'`
				if [ ! -z "$prot" ]; then
	                zcat $DATABASE"/proteins/Macaca_mulatta.faa.gz" | awk -v RS="(^|\n)>" '$1=="'$prot'" {print ">"$0}' >>  $PROJECT"/results/seqs/"$NP_number"_ClinVar.fasta"
	                echo $NP_number "-->" $gene_id "-->" $id "-->" $prot "[for Macaque]" >> "logs/get_homology_seq_clinvar.log"
				fi
	    	done
			#Rat
			echo $homology_line | awk -v FS=";" '{print $4}' | tr ',' '\n' | while read id
	        do
	            prot=`cat $DATABASE"/mapping/NP_tables/Rattus_norvegicus.list" | awk -v FS="\t" '$6=="'$id'" {print $6,$8,$9}' | sort -rnk3,3 | awk 'NR==1 {print $2}'`
	            if [ ! -z "$prot" ]; then
					zcat $DATABASE"/proteins/Rattus_norvegicus.faa.gz" | awk -v RS="(^|\n)>" '$1=="'$prot'" {print ">"$0}' >>  $PROJECT"/results/seqs/"$NP_number"_ClinVar.fasta"
	                echo $NP_number "-->" $gene_id "-->" $id "-->" $prot "[for Rat]" >> "logs/get_homology_seq_clinvar.log"
				fi
	        done
			#Mouse
			echo $homology_line | awk -v FS=";" '{print $5}' | tr ',' '\n' | while read id
	        do
	            prot=`cat $DATABASE"/mapping/NP_tables/Mus_musculus.list" | awk -v FS="\t" '$6=="'$id'" {print $6,$8,$9}' | sort -rnk3,3 | awk 'NR==1 {print $2}'`
	            if [ ! -z "$prot" ]; then
					zcat $DATABASE"/proteins/Mus_musculus.faa.gz" | awk -v RS="(^|\n)>" '$1=="'$prot'" {print ">"$0}' >>  $PROJECT"/results/seqs/"$NP_number"_ClinVar.fasta"
	                echo $NP_number "-->" $gene_id "-->" $id "-->" $prot "[for Mouse]" >> "logs/get_homology_seq_clinvar.log"
				fi
	        done
			#ZebraFish
			echo $homology_line | awk -v FS=";" '{print $6}' | tr ',' '\n' | while read id
	        do
	        	prot=`cat $DATABASE"/mapping/NP_tables/Danio_rerio.list" | awk -v FS="\t" '$6=="'$id'" {print $6,$8,$9}' | sort -rnk3,3 | awk 'NR==1 {print $2}'`
	        	if [ ! -z "$prot" ]; then
					zcat $DATABASE"/proteins/Danio_rerio.faa.gz" | awk -v RS="(^|\n)>" '$1=="'$prot'" {print ">"$0}' >>  $PROJECT"/results/seqs/"$NP_number"_ClinVar.fasta"
	                echo $NP_number "-->" $gene_id "-->" $id "-->" $prot "[for ZebraFish]" >> "logs/get_homology_seq_clinvar.log"
				fi
	        done
			#Frog
			echo $homology_line | awk -v FS=";" '{print $7}' | tr ',' '\n' | while read id
	        do
	            prot=`cat $DATABASE"/mapping/NP_tables/Xenopus_tropicalis.list" | awk -v FS="\t" '$6=="'$id'" {print $6,$8,$9}' | sort -rnk3,3 | awk 'NR==1 {print $2}'`
				if [ ! -z "$prot" ]; then
	            	zcat $DATABASE"/proteins/Xenopus_tropicalis.faa.gz" | awk -v RS="(^|\n)>" '$1=="'$prot'" {print ">"$0}' >>  $PROJECT"/results/seqs/"$NP_number"_ClinVar.fasta"
	                echo $NP_number "-->" $gene_id "-->" $id "-->" $prot "[for Frog]" >> "logs/get_homology_seq_clinvar.log"
				fi
	        done
			#Fruitfly
			echo $homology_line | awk -v FS=";" '{print $8}' | tr ',' '\n' | while read id
	        do
	            prot=`cat $DATABASE"/mapping/NP_tables/Drosophila_melanogaster.list" | awk -v FS="\t" '$6=="'$id'" {print $6,$8,$9}' | sort -rnk3,3 | awk 'NR==1 {print $2}'`
				if [ ! -z "$prot" ]; then
	            	zcat $DATABASE"/proteins/Drosophila_melanogaster.faa.gz" | awk -v RS="(^|\n)>" '$1=="'$prot'" {print ">"$0}' >>  $PROJECT"/results/seqs/"$NP_number"_ClinVar.fasta"
	                echo $NP_number "-->" $gene_id "-->" $id "-->" $prot "[for Fruitfly]" >> "logs/get_homology_seq_clinvar.log"
				fi
	        done
			#Celegans
			echo $homology_line | awk -v FS=";" '{print $9}' | tr ',' '\n' | while read id
	        do
	            prot=`cat $DATABASE"/mapping/NP_tables/Caenorhabditis_elegans.list" | awk -v FS="\t" '$6=="'$id'" {print $6,$8,$9}' | sort -rnk3,3 | awk 'NR==1 {print $2}'`
				if [ ! -z "$prot" ]; then
	            	zcat $DATABASE"/proteins/Caenorhabditis_elegans.faa.gz" | awk -v RS="(^|\n)>" '$1=="'$prot'" {print ">"$0}' >>  $PROJECT"/results/seqs/"$NP_number"_ClinVar.fasta"
	                echo $NP_number "-->" $gene_id "-->" $id "-->" $prot "[for Celegans]" >> "logs/get_homology_seq_clinvar.log"
				fi
	        done
		fi
	done
}

#GET HOMOLOGY SEQ FOR gnomAD
get_homology_seq_gnomAD(){
	if [ -d "$PROJECT/results/seqs2" ]; then rm -r $PROJECT"/results/seqs2"; fi
	cp $PROJECT"/results/gnomad_seqs" $PROJECT"/results/seqs2" -r 
	if [ -f "logs/get_homology_seq_gnomad.log" ]; then rm "logs/get_homology_seq_gnomad.log"; fi

	ls -l /opt/current_project/results/gnomad_seqs/*_gnomAD.fasta | awk -F"/" '{print $6}' | awk -F"_gnomAD.fasta" '{print $1}' > $PROJECT"/tmp/tmp_gnomAD.list"
	IFS=$'\n'; for ENST_number in $(cat "$PROJECT/tmp/tmp_gnomAD.list"); 
	do
		gene_id=`cat $DATABASE"/mapping/NewCurated_ENSTvsGENEID.csv" | awk -v FS="," '$1=="'$ENST_number'" {print $2}' | sort -u | head -n 1`
		#HOMOLOGY PART STARTS HERE
		if [ -z "$gene_id" ]; then
			echo "GeneID is empty! for $ENST_number" >> "logs/get_homology_seq_gnomad.log";
		else
			homology_line=`cat $DATABASE"/orthology/Homology.list" | awk -v FS="\t" -v OFS=";" '$1=="'$gene_id'" {print $1,$2,$3,$4,$5,$6,$7,$8,$9}'`
			#Chimp
			echo $homology_line | awk -v FS=";" '{print $2}' | tr ',' '\n' | while read id
			do
				prot=`cat $DATABASE"/mapping/NP_tables/Pan_troglodytes.list" | awk -v FS="\t" '$6=="'$id'" {print $6,$8,$9}' | sort -rnk3,3 | awk 'NR==1 {print $2}'`
				if [ ! -z "$prot" ]; then
					zcat $DATABASE"/proteins/Pan_troglodytes.faa.gz" | awk -v RS="(^|\n)>" '$1=="'$prot'" {print ">"$0}' >>  $PROJECT"/results/seqs2/"$ENST_number"_gnomAD.fasta"
					echo $ENST_number "-->" $gene_id "-->" $id "-->" $prot "[for Chimp]" >> "logs/get_homology_seq_gnomad.log"
				fi
			done
			#Macaque
			echo $homology_line | awk -v FS=";" '{print $3}' | tr ',' '\n' | while read id
	        do
	            prot=`cat $DATABASE"/mapping/NP_tables/Macaca_mulatta.list" | awk -v FS="\t" '$6=="'$id'" {print $6,$8,$9}' | sort -rnk3,3 | awk 'NR==1 {print $2}'`
				if [ ! -z "$prot" ]; then
	                zcat $DATABASE"/proteins/Macaca_mulatta.faa.gz" | awk -v RS="(^|\n)>" '$1=="'$prot'" {print ">"$0}' >>  $PROJECT"/results/seqs2/"$ENST_number"_gnomAD.fasta"
	                echo $ENST_number "-->" $gene_id "-->" $id "-->" $prot "[for Macaque]" >> "logs/get_homology_seq_gnomad.log"
				fi
	    	done
			#Rat
			echo $homology_line | awk -v FS=";" '{print $4}' | tr ',' '\n' | while read id
	        do
	            prot=`cat $DATABASE"/mapping/NP_tables/Rattus_norvegicus.list" | awk -v FS="\t" '$6=="'$id'" {print $6,$8,$9}' | sort -rnk3,3 | awk 'NR==1 {print $2}'`
	            if [ ! -z "$prot" ]; then
					zcat $DATABASE"/proteins/Rattus_norvegicus.faa.gz" | awk -v RS="(^|\n)>" '$1=="'$prot'" {print ">"$0}' >>  $PROJECT"/results/seqs2/"$ENST_number"_gnomAD.fasta"
	                echo $ENST_number "-->" $gene_id "-->" $id "-->" $prot "[for Rat]" >> "logs/get_homology_seq_gnomad.log"
				fi
	        done
			#Mouse
			echo $homology_line | awk -v FS=";" '{print $5}' | tr ',' '\n' | while read id
	        do
	            prot=`cat $DATABASE"/mapping/NP_tables/Mus_musculus.list" | awk -v FS="\t" '$6=="'$id'" {print $6,$8,$9}' | sort -rnk3,3 | awk 'NR==1 {print $2}'`
	            if [ ! -z "$prot" ]; then
					zcat $DATABASE"/proteins/Mus_musculus.faa.gz" | awk -v RS="(^|\n)>" '$1=="'$prot'" {print ">"$0}' >>  $PROJECT"/results/seqs2/"$ENST_number"_gnomAD.fasta"
	                echo $ENST_number "-->" $gene_id "-->" $id "-->" $prot "[for Mouse]" >> "logs/get_homology_seq_gnomad.log"
				fi
	        done
			#ZebraFish
			echo $homology_line | awk -v FS=";" '{print $6}' | tr ',' '\n' | while read id
	        do
	        	prot=`cat $DATABASE"/mapping/NP_tables/Danio_rerio.list" | awk -v FS="\t" '$6=="'$id'" {print $6,$8,$9}' | sort -rnk3,3 | awk 'NR==1 {print $2}'`
	        	if [ ! -z "$prot" ]; then
					zcat $DATABASE"/proteins/Danio_rerio.faa.gz" | awk -v RS="(^|\n)>" '$1=="'$prot'" {print ">"$0}' >>  $PROJECT"/results/seqs2/"$ENST_number"_gnomAD.fasta"
	                echo $ENST_number "-->" $gene_id "-->" $id "-->" $prot "[for ZebraFish]" >> "logs/get_homology_seq_gnomad.log"
				fi
	        done
			#Frog
			echo $homology_line | awk -v FS=";" '{print $7}' | tr ',' '\n' | while read id
	        do
	            prot=`cat $DATABASE"/mapping/NP_tables/Xenopus_tropicalis.list" | awk -v FS="\t" '$6=="'$id'" {print $6,$8,$9}' | sort -rnk3,3 | awk 'NR==1 {print $2}'`
				if [ ! -z "$prot" ]; then
	            	zcat $DATABASE"/proteins/Xenopus_tropicalis.faa.gz" | awk -v RS="(^|\n)>" '$1=="'$prot'" {print ">"$0}' >>  $PROJECT"/results/seqs2/"$ENST_number"_gnomAD.fasta"
	                echo $ENST_number "-->" $gene_id "-->" $id "-->" $prot "[for Frog]" >> "logs/get_homology_seq_gnomad.log"
				fi
	        done
			#Fruitfly
			echo $homology_line | awk -v FS=";" '{print $8}' | tr ',' '\n' | while read id
	        do
	            prot=`cat $DATABASE"/mapping/NP_tables/Drosophila_melanogaster.list" | awk -v FS="\t" '$6=="'$id'" {print $6,$8,$9}' | sort -rnk3,3 | awk 'NR==1 {print $2}'`
				if [ ! -z "$prot" ]; then
	            	zcat $DATABASE"/proteins/Drosophila_melanogaster.faa.gz" | awk -v RS="(^|\n)>" '$1=="'$prot'" {print ">"$0}' >>  $PROJECT"/results/seqs2/"$ENST_number"_gnomAD.fasta"
	                echo $ENST_number "-->" $gene_id "-->" $id "-->" $prot "[for Fruitfly]" >> "logs/get_homology_seq_gnomad.log"
				fi
	        done
			#Celegans
			echo $homology_line | awk -v FS=";" '{print $9}' | tr ',' '\n' | while read id
	        do
	            prot=`cat $DATABASE"/mapping/NP_tables/Caenorhabditis_elegans.list" | awk -v FS="\t" '$6=="'$id'" {print $6,$8,$9}' | sort -rnk3,3 | awk 'NR==1 {print $2}'`
				if [ ! -z "$prot" ]; then
	            	zcat $DATABASE"/proteins/Caenorhabditis_elegans.faa.gz" | awk -v RS="(^|\n)>" '$1=="'$prot'" {print ">"$0}' >>  $PROJECT"/results/seqs2/"$ENST_number"_gnomAD.fasta"
	                echo $ENST_number "-->" $gene_id "-->" $id "-->" $prot "[for Celegans]" >> "logs/get_homology_seq_gnomad.log"
				fi
	        done
		fi
	done

	rm $PROJECT"/tmp/tmp_gnomAD.list"
}
get_homology_seq_gnomAD

#NEW CONVERSION TABLE
conversion(){
	IFS=$'\n'; for ENST_number in $(cat "$DATABASE/mapping/gnomAD_ENST_unique.list"); 
	do
		hgnc_id=`cat $DATABASE"/mapping/ENST_ENSTV_GeneID.list" | awk -F"," '$1=="'$ENST_number'" {print $4}' | sort -u`
		hgnc_id_number=`cat $DATABASE"/mapping/ENST_ENSTV_GeneID.list" | awk -F"," '$1=="'$ENST_number'" {print $4}' | sort -u | wc -l`
		if [ ! -z "$hgnc_id" ] && [ "$hgnc_id_number" == "1" ]; then
			term="HGNC:$hgnc_id"
			gene_id=`cat $DATABASE"/mapping/HGNC.list" | awk -F"\t" '$1=="'$term'" {print $5}' | sort -u | head -n 1`
		else
			missing_gene_name=`cat $DATABASE"/mapping/ENST_ENSTV_GeneID.list" | awk -F"," '$1=="'$ENST_number'" {print $3}' | sort -u`
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

<< idk
echo "Clinvar is running now"
get_homology_seq_clinvar &
echo "gnomAD is running now"
get_homology_seq_gnomAD &
wait $(jobs -p)
python3 create_alignments.py &
python3 db_mapping.py &
wait $(jobs -p)
python3 compute_conservation_scores.py
idk