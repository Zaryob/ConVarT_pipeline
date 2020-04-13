BEGIN{
	RS="//"
	FS="\n"
	nid=0;
	tempversion="NA";
	say=1;
	synend="|"
}
{
	
	for (i=1;i<=NF;i++)
{
	
	locps=index($i,"LOCUS")+5;
	accps=index($i,"ACCESSION")+9;
	verps=index($i,"VERSION")+7;

	dbps=index($i,"/db_xref=")+9;
	geneps=index($i,"/gene=")+6;
	genesynps=index($i,"/gene_synonym=")+14;
	stdnameps=index($i,"/standard_name=")+15;
	locustagps=index($i,"/locus_tag=")+11;
	
	#print 	locps,accps,verps,dbps,geneps,genesynps;
		if(locps>5)
		{
				#print ">>" NR,$i;
				split($i,alocus," ");
				seqlen=alocus[3];
		}

	gsub(/^[ \t]+/,"",$i);
	
	gsub(" ","",$i);
	veri=substr($i,1);
	$i=substr($i,2);



	if(accps>9)
		{
				#print ">>" NR,$i;
		}

	if(verps>7)
		{
				versionid=substr(veri,verps);
				#print $i,versionid
				gsub(/^[ \t]+/,"",versionid);

				if (tempversion!=versionid)
				{
					say=1;
					if (NR>1)
					{

						printf tempversion" db_xrefs=";
						for(di=1;di<length(db_xref[NR-1]);di++)
							{

								printf substr(db_xref[NR-1][di],9)",";
							}

							printf  substr(db_xref[NR-1][di],9) "\n";

					}

					tempversion=versionid;

				}

				print versionid,"sequence_length="seqlen;



		}

	if(dbps>9)
		{			
				
				db_xref[NR][say]=$i;
				say++;
		}

	if(geneps>6)
		{
				print versionid,$i;
		}

	if(genesynps>14 || synend==";")
		{
				if(synend!=";")
{
				printf versionid" "$i
}

				if(synend==";")
{
                                printf $i
}                               

		
	synend=substr($i,length($i))
	if (synend!=";")
		{    
			printf "\n"
			}
		#print ">>>>>>>",synend		
		}

	if(stdnameps>15)
		{
				print versionid,$i;
		}
	if(locustagps>11)
		{
				print versionid,$i;
		}

	locps=0;
	accps=0;
	verps=0;

	dbps=0;
	geneps=0;
	genesynps=0;
	stdnameps=0;
	locustagps=0;


}	

}
