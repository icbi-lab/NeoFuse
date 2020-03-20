#! /bin/bash
OUTDIR=$PWD
declare -i PEPMAX
declare -i PEPMIN
PEPMIN=8
declare -i CORES
CORES=1
THRESHOLD=""
RANK=""
CONF="L"
NETMHCPAN="false"
BUILD="false"
BUILDREF="false"
VERSION="GRCh38"

while getopts ":1::2::d::i::o::m::M::n::t::T::c::s::g::a::h" opt;
do
	case $opt in
	1)	READ1="$OPTARG";;
	2)	READ2="$OPTARG";;
	d)	ID="$OPTARG";;
	i)	IN="$OPTARG";;
	o)  OUTDIR="$OPTARG";;
	m)	PEPMIN="$OPTARG";;
	M)	PEPMAX="$OPTARG";;
	n)	CORES="$OPTARG";;
	t)	THRESHOLD="$OPTARG";;
	T)	RANK="$OPTARG";;
	c)	CONF="$OPTARG";;
	s)	STARINDEX="$OPTARG";;
	g)	GENOMEDIR="$OPTARG";;
	a)	ANNOTATION="$OPTARG";;
	N)	NETMHCPAN="$OPTARG";;
	h)	HELP="true";;
	esac
done

if [ "$READ1" == "" ] && [ "$IN" == "" ] && [ "$BUILDREF" == "false" ] || [ "$HELP" == "true" ]; then
	echo "Usage: NeoFuse <arguments> [options]

	<Arguments>

	-1: Path to read 1 FASTQ file (mandatory only for single sample analysis)
	-2: Path to read 2 FASTQ fie (optional for single-end reads)
	-i: Path to TSV file (mandatory only for multiple sample analysis)
	-s: Path to STAR index directory
	-g: Path to reference genome FASTA file 
	-a: Path to annotation GTF file
	-o: Output directory (default: "./")
	-d: Run ID (default: input filename)

	-B: Build singularity/docker image
	-R: Build indexes

	[Options]

	-m: Minimum peptide length (8, 9, 10, or 11; default: 8)
	-M: Maximum peptide length (8, 9, 10, or 11; default: 8)
	-n: Number of cores (default: 1)
	-t: IC50 binding affinity threshold (default: 500)
	-T: Percentile rank threshold (default: Inf)
	-c: Mimimum confidence score (H, M or L; default: L)

	--singularity: NeoFuse will use the Singularity image
	--docker: NeoFuse will use the Docker image
	
	For more information please refer to https://icbi.i-med.ac.at/software/NeoFuse/
	"
	exit 0
else
	:
fi

if [ "$PEPMAX" == "" ]; then
	PEPMAX=$PEPMIN
fi

if [ "$IN" != "" ] && [ "$READ1" != "" ]; then
	echo " You may specify as input a TSV file OR FASTQ file(s), but not both "
	echo "Exiting ..."
	exit 1
fi

if [ "$IN" != "" ] && [ "$READ2" != "" ]; then
	echo " You may specify as input a TSV file OR FASTQ file(s), but not both "
	echo "Exiting ..."
	exit 1
fi

if [ "$STARINDEX" == "" ]; then
	echo "You must specify the STAR index directory"
	echo "Exiting ..."
	exit 1
elif [ "$GENOMEDIR" == "" ]; then
	echo "You must specify a FASTA reference file"
	echo "Exiting ..."
	exit 1
elif [ "$ANNOTATION" == "" ]; then
	echo "You must specify an annotation GTF file"
	echo "Exiting ..."
	exit 1
fi

OUTDIR=`realpath $OUTDIR`

# Check for multiple or single sample(s) - call the appropriate script
if [ "$IN" != "" ]; then
	NeoFuse_multi -i $IN -o $OUTDIR -m $PEPMIN -M $PEPMAX -n $CORES -t $THRESHOLD -T $RANK -c $CONF -s $STARINDEX -g $GENOMEDIR -a $ANNOTATION -N $NETMHCPAN -r $OUTDIR
elif [ "$READ1" != "" ]; then
	NeoFuse_single -1 $READ1 -2 $READ2 -d $ID -o $OUTDIR -m $PEPMIN -M $PEPMAX -n $CORES -t $THRESHOLD -T $RANK -c $CONF -s $STARINDEX -g $GENOMEDIR -a $ANNOTATION -N $NETMHCPAN -r $OUTDIR
else
	echo "Invalid arguments"
	echo "Please check your input files"
	echo "Exiting ..."
	exit 1
fi
