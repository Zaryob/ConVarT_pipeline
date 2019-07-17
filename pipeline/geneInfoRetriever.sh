#!/bin/bash
PROJECT="/opt/current_project"
DATABASE="/opt/current_project/db"
#Human_GPFF_Formatted.list

echo "prot_id,seq_len,gene_symbol,genesynonyms,tax_id,gene_id,HGNC_id,MIM_id"
#cat $DATABASE"/others/bilinen_genler.fasta" | awk '/>/ {print $2}' | while read np_number
cat $DATABASE"/others/Human_GPFF_Formatted_ProtIDs.txt" | while read np_number
do
        line=$(cat $DATABASE"/others/Human_GPFF_Formatted.list" | awk '$1=="'$np_number'" {print $2}' | tr '\n' '|' | awk '{print "\n"$0'})
        slen=$(echo $line | awk -F"[|=]" '{print $2}')
        gsymbol=$(echo $line | awk -F"[|=]" '{print $4}' | tr -d '"')
        gsynonyms=$(echo $line | awk -F"[|=]" '{print $6}' | tr -d '"')
        taxid=$(echo $line | awk -F"[|=]" '{print $8}' | tr -d '"' | awk -F'[:,]' '{print $2}')
        geneid=$(echo $line | awk -F"[|=]" '{print $8}' | tr -d '"' | awk -F'[:,]' '{print $4}')
        hgncid=$(echo $line | awk -F"[|=]" '{print $8}' | tr -d '"' | awk -F'[:,]' '{print $7}')
        mimid=$(echo $line | awk -F"[|=]" '{print $8}' | tr -d '"' | awk -F'[:,]' '{print $9}')
	oinfo=$(echo $line | awk -F"[|=]" '{print $8}' | tr -d '"' | tr ',' ';')

        seq_len=$([ -z "$slen" ] && echo "Empty" || echo $slen)
        gene_symbol=$([ -z "$gsymbol" ] && echo "Empty" || echo $gsymbol)
        gene_synonyms=$([ -z "$gsynonyms" ] && echo "Empty" || echo $gsynonyms)
	other_info=$([ -z "$oinfo" ] && echo "Empty" || echo $oinfo)

	echo $np_number","$seq_len","$gene_symbol","$gene_synonyms","$other_info
#        tax_id=$([ -z "$taxid" ] && echo "Empty" || echo $taxid)
#        gene_id=$([ -z "$geneid" ] && echo "Empty" || echo $geneid)
#        hgnc_id=$([ -z "$hgncid" ] && echo "Empty" || echo $hgncid)
#        mim_id=$([ -z "$mimid" ] && echo "Empty" || echo $mimid)
#        echo $np_number","$seq_len","$gene_symbol","$gene_synonyms","$tax_id","$gene_id","$hgnc_id","$mim_id
done

