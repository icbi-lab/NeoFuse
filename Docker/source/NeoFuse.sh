#! /bin/bash
OUTDIR="./"
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

while getopts ":1::2::d::i::o::m::M::n::t::T::c::s::g::a::" opt;
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
	esac
done

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

# Check for multiple or single sample(s) - call the appropriate script
if [ "$IN" != "" ]; then
	NeoFuse_multi -i $IN -o $OUTDIR -m $PEPMIN -M $PEPMAX -n $CORES -t $THRESHOLD -T $RANK -c $CONF -s $STARINDEX -g $GENOMEDIR -a $ANNOTATION -N $NETMHCPAN
elif [ "$READ1" != "" ]; then
	NeoFuse_single -1 $READ1 -2 $READ2 -d $ID -o $OUTDIR -m $PEPMIN -M $PEPMAX -n $CORES -t $THRESHOLD -T $RANK -c $CONF -s $STARINDEX -g $GENOMEDIR -a $ANNOTATION -N $NETMHCPAN
else
	echo "Invalid arguments"
	echo "Please check your input files"
	echo "Exiting ..."
	exit 1
fi