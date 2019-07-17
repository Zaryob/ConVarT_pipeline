curdir=$(pwd)
targetdir=$curdir/$1
inputdir=$curdir/$2
echo $inputdir
ls $inputdir > $curdir/"inputfiles.txt"

cat $curdir/inputfiles.txt | while read inputfile
	do
		echo $inputfile
		echo $inputdir"/"$inputfile
		filestring=$(basename $inputfile)
		organism=${filestring:0:-8}

echo $targetdir/$organism"."$2".list"
echo $inputdir/$inputfile

		awk -f /mnt/depo/projects/ciliogenics/scripts/gpff_ret.awk \
		<(zcat $inputdir"/"$inputfile) \
		> $targetdir"/"$organism"."$2".list"

	done
