
for i in *_unassigned.txt
do
	# Get basename
	basename=`echo $i | cut -d _ -f1,2,3,4`
	temp_fastq_name="original_fastqs/""$basename""_001.fastq.gz"
	echo $temp_fastq_name
	# Pull readnames
	out="$basename.readnames.txt"
	cat $i | grep ">" | tr -d ">" > $out

	# Pull reads from orignal fastqs 
	 echo "retrieving original fastq records for $i"
	 in="$temp_fastq_name"
	 out="$basename.unassigned.fastq"
	seqtk subseq $in $basename.readnames.txt > $out

	# Retrim files
	 echo "trimming $i"
	 in="$basename.unassigned.fastq"
	 out="$basename.trimmed.fastq"
	 trimmomatic SE -phred33 $in $out -threads 8 ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

	# #remove mitochondrial sequence

	 echo "depleting mitochondria for $i"
	 in="$basename.trimmed.fastq "
	 out="$basename.mito_removed.fastq"
	 bowtie2 -x bt2/mito -p 8 -U $in | samtools view -Sb -f 4 - | samtools fastq - > $out

	# make pseudo paired fastq
	mkdir $basename
	 in="in=""$basename.mito_removed.fastq"
	 rone="out1=""$basename.r1.fastq"
	 rtwo="out2=""$basename.r2.fastq"

	 bbfakereads.sh $in $rone $rtwo
	 mv $basename.r1.fastq $basename/
	 mv $basename.r2.fastq $basename/
done


# spades command 
# spades -1 *r1* -2 *r2* -t 42  --meta -k 125,127 -o <assembled_name> 
