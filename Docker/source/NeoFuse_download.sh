#! /bin/bash

CORES=1
VERSION="GRCh38"
while getopts "o::t::v::" opt; 
do
	case $opt in
	o)  OUT="$OPTARG";;
	t)	THREADS="$OPTARG";;
	v)	VERSION="$OPTARG";;
	esac
done

#OUT=`realpath $OUT`

INDEXDIR=$OUT"/STAR_idx"

mkdir -p $INDEXDIR

cd $OUT

if test -z "$OUT" 
then
	echo "Usage: genomes_download -o [output] -t [cores] -v [genome version]"
	echo "You must specify an output directory"
	echo "Exiting ..."
	exit 1
fi

if [ $VERSION == "GRCh38" ]; then
	wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.annotation.gtf.gz
	gunzip gencode.v31.annotation.gtf.gz
	wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/GRCh38.primary_assembly.genome.fa.gz
	gunzip GRCh38.primary_assembly.genome.fa.gz

	STAR --runMode genomeGenerate --runThreadN $THREADS --genomeDir $INDEXDIR \
	--genomeFastaFiles $OUT"/GRCh38.primary_assembly.genome.fa" \
	--sjdbGTFfile $OUT"/gencode.v31.annotation.gtf"
	if [ `echo $?` != 0 ]; then
		exit 1
	else
		:
	fi
	echo "|---------------------------------------------------------"
	echo "| genomes_download has succesfully finished			   "
	echo "| Genome Version: GRCh38.p12 								"
	echo "| FASTA File: $OUT\"/GRCh38.primary_assembly.genome.fa\"  "
	echo "| GTF File: $OUT\"/gencode.v31.annotation.gtf\"		   "
	echo "| STAR Index: $INDEXDIR									"
	echo "|---------------------------------------------------------"
elif [ $VERSION == "GRCh37" ]; then
	wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
	gunzip gencode.v19.annotation.gtf.gz
	wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz
	gunzip GRCh37.p13.genome.fa.gz

	STAR --runMode genomeGenerate --runThreadN $THREADS --genomeDir $INDEXDIR \
	--genomeFastaFiles $OUT"/GRCh37.p13.genome.fa" \
	--sjdbGTFfile $OUT"/gencode.v19.annotation.gtf"
	if [ `echo $?` != 0 ]; then
		exit 1
	else
		:
	fi
	echo "|--------------------------------------------------------"
	echo "| genomes_download has succesfully finished			  "
	echo "| Genome Version: GRCh37.p13 							   "
	echo "| FASTA File: $OUT\"/GRCh37.p13.genome.fa\"			  "
	echo "| GTF File: $OUT\"/gencode.v19.annotation.gtf\"		  "
	echo "| STAR Index: $INDEXDIR								   "
	echo "|--------------------------------------------------------"
else
	echo "|--------------------------------|"
	echo "| Error: invalid genome version  |"
	echo "| Please choose GRCh37 or GRCh38 |"
	echo "|--------------------------------|"
	exit 1
fi
