https://www.unix.com/shell-programming-and-scripting/217359-help-awk-sed-putting-space-after-numbers-separate-number-characters.html
awk -v FS="\t" '{print $8}' firstFormat2Şub2019.csv | awk -F"[().]" '{print $1,$3,$7}' |head -n 3
 sed 's/^[0-9]*/& /'  file
ÇALŞAN KOMUT:::: echo "abce1234jklm" | sed 's/[0-9][0-9]*/ & /'
https://www.biostars.org/p/137256/
#kotu https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id=4&retmode=xml
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&id=4&rettype=variation
