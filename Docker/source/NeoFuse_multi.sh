#! /bin/bash

BLACKLIST=/usr/local/bin/genomes/blacklist_hg38_GRCh38_2018-11-04.tsv.gz
OPTIREF=/usr/local/bin/OptiType-1.3.2/data/hla_reference_rna.fasta
Optitype=/usr/local/bin/OptiType-1.3.2/OptiTypePipeline.py
YARAIDX=/usr/local/bin/yara_idx/hla_reference_rna
OUTDIR=$PWD
declare -i PEPMAX
declare -i PEPMIN
PEPMIN=8
declare -i CORES
CORES=1
declare -i RAMLIMIT
RAMLIMIT=0
THRESHOLD=""
RANK=""
CONF="L"
CUSTOMLIST="false"
KEEPBAM="false"
NETMHCPAN="false"

while getopts "i:o::m::M::n::t::T::c::s::g::a::r::C::N::l::k::" opt;
do
	case $opt in
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
	r)	REALOUT="$OPTARG";;
	l)	RAMLIMIT="$OPTARG";;
	C)	CUSTOMLIST="$OPTARG";;
	k)	KEEPBAM="$OPTARG";;
	N)	NETMHCPAN="$OPTARG";;
	esac
done

# check if MHCflurry models exist
if [ ! -d "~/.local/share/mhcflurry/4/2.0.0/models_class1/" ]; then
	mhcflurry-downloads fetch &> /dev/null
fi

# Check input TSV file
echo > $OUTDIR"/in.tsv"
while read -r line
do
	if [[ $line =~ \# ]]; then
		:
	elif [[ $line =~ " " ]]; then
		echo "No whitespaces allowed in the input TSV file"
		rm $OUTDIR"/in.tsv"
		echo "Exiting ..."
		exit 1
	else
		echo $line >> $OUTDIR"/in.tsv"
	fi
done < $IN

# process args

if test -z "$REALOUT"
then
	REALOUT=$OUTDIR
fi

if test -z "$IN"
then
	echo "usage: NeoFuse_multi -i <input> -o [output] -m [peptide min length] -M [peptide Max length] -n [cores] -c [confidence level] -t [affinity threshold] -s <STAR index> -g <FASTA file> -a <GTF file>"
	echo "You must specify an input file"
	echo "Exiting ..."
	exit 1
fi

if [ "$THRESHOLD" == "" ] && [ "$RANK" == "" ]; then
	THRESHOLD=500
	RANK=inf
elif [ "$THRESHOLD" != "" ] && [ "$RANK" == "" ]; then
	RANK=inf
elif [ "$THRESHOLD" == "" ] && [ "$RANK" != "" ]; then
	THRESHOLD=inf
fi

if test -z "$PEPMAX"
then
	PEPMAX=$PEPMIN
fi

if [ $PEPMIN -lt 8 ] || [ $PEPMIN -gt 11 ]; then
	echo "Invalid option -m [must be 8,9,10 or 11]"
	exit 1
fi

if [ $PEPMAX -lt 8 ] || [ $PEPMAX -gt 11 ]; then
	echo "Invalid option -M [must be 8,9,10 or 11]"
	exit 1
fi

if [ $PEPMIN -eq $PEPMAX ]; then
	:
elif [ $PEPMIN -gt $PEPMAX ]; then
	echo "Invalid options -m and -M [ m > M ]"
	exit 1
fi

PEPLEN=( )
for i in `seq $PEPMIN $PEPMAX`; do
	PEPLEN+="$i "
done

case "$CONF" in
	"H") CONFIDENCE="high" ;;
	"M") CONFIDENCE="high medium" ;;
	"L") CONFIDENCE="high medium low" ;;
esac

# Manage cores and multithreading
if [ $CORES -le 3 ]; then
	STARTHREADS=1
	ARRIBATHREADS=1
	RAZERTHREADS=1
else
	if [ $(($CORES % 2)) -eq 0 ]; then
		STARTHREADS=$(($CORES / 2))
		ARRIBATHREADS=$(($CORES / 2))
		RAZERTHREADS=$(($CORES / 2))
	else
		STARTHREADS=$((($CORES -1) / 2))
		ARRIBATHREADS=$((($CORES -1) / 2))
		RAZERTHREADS=$((($CORES -1) / 2))
	fi
fi

# Read the TSV file, skip 1st line, get vars
## Process each pair of reads
sed 1d $OUTDIR"/in.tsv" | while read -r id read1 read2
do
	FILE=$id
	READ1=$read1
	READ2=$read2
	# Check for user defined ID, if not get ID from file names
	if test $FILE; then
		:
	else
		FILE=$(echo ${READ1##*/} | sed s/_.*fastq.*//)
	fi

	# Get file names - to be used for processing
	NAME1=$(echo ${READ1##*/})
	NAME2=$(echo ${READ2##*/})

	# Check if FASTQ files are zipped - if yes issue zcat argument to STAR
	if [[ $NAME1 == *".gz"* ]]; then
	  ZIPPED="--readFilesCommand zcat"
	fi

	# Create output folders
	OUTDIRALIGN=$OUTDIR"/"$FILE"/STAR/"
	OUTDIRARRIBA=$OUTDIR"/"$FILE"/Arriba/"
	OUTDIRCOUNTS=$OUTDIR"/"$FILE"/FeatureCounts/"
	OUTDIRTPM=$OUTDIR"/"$FILE"/TPM/"
	OUTDIRRPKM=$OUTDIR"/"$FILE"/RPKM/"
	TEMPDIROPTI=$OUTDIR"/"$FILE"/OptiType/tmp/"
	OUTDIROPTI=$OUTDIR"/"$FILE"/OptiType/"
	OUTDIRCLEAVEPEP=$OUTDIR"/"$FILE"/Peptides/"
	FINALOUTDIR=$OUTDIR"/"$FILE"/NeoFuse/"
	FINALTMP=$OUTDIR"/"$FILE"/NeoFuse/tmp/"
	LOGSDIR=$OUTDIR"/"$FILE"/LOGS/"
	mkdir -p $OUTDIRALIGN
	mkdir -p $OUTDIRARRIBA
	mkdir -p $OUTDIRCOUNTS
	mkdir -p $OUTDIRTPM
	mkdir -p $OUTDIRRPKM
	mkdir -p $OUTDIROPTI
	mkdir -p $TEMPDIROPTI
	mkdir -p $OUTDIRCLEAVEPEP
	mkdir -p $FINALOUTDIR
	mkdir -p $FINALTMP
	mkdir -p $LOGSDIR

	# Check for PE or SE reads and process the files accordingly
	if test -f "$READ2"; then
		echo "[-------------------------------- [NeoFuse] --------------------------------]"
		echo
		echo " Paired End (PE) Reads detected: commencing processing" | sed "s/^/[NeoFuse] /"
		echo " Processing files $NAME1 - $NAME2" | sed "s/^/[NeoFuse] /"
		echo " STAR Run started at:" `date +"%T"` | sed "s/^/[NeoFuse] /"
		# STAR
		STAR --runThreadN $STARTHREADS \
		--runMode alignReads \
		--genomeDir $STARINDEX \
		--readFilesIn $READ1 $READ2 $ZIPPED \
		--outFileNamePrefix $OUTDIRALIGN$FILE"." \
		--genomeLoad NoSharedMemory \
		--outReadsUnmapped Fastx \
		--outStd BAM_Unsorted \
		--outBAMcompression 0 \
		--outSAMtype BAM Unsorted 2>$LOGSDIR$FILE.STAR.err |
			samtools sort -@ $SAMTOOLSTHREADS -O BAM -o ${OUTDIRALIGN}${FILE}.Aligned.sortedByCoord.out.bam /dev/stdin \
			> $LOGSDIR$FILE.samtools.log 2>$LOGSDIR$FILE.samtools.err &
		

		# Arriba
		echo " Arriba Run started at:" `date +"%T"` | sed "s/^/[NeoFuse] /"
		STAR --runThreadN $ARRIBATHREADS \
		--genomeDir $STARINDEX --genomeLoad NoSharedMemory --limitBAMsortRAM $RAMLIMIT \
		--readFilesIn $READ1 $READ2 $ZIPPED \
		--outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 \
		--outFilterMultimapNmax 1 --outFilterMismatchNmax 3 \
		--outFileNamePrefix $OUTDIRARRIBA$FILE"." \
		--chimSegmentMin 10 --chimOutType WithinBAM SoftClip --chimJunctionOverhangMin 10 --chimScoreMin 1 \
		--chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --alignSJstitchMismatchNmax 5 -1 5 5 \
		--chimSegmentReadGapMax 3 |
		arriba \
		-x /dev/stdin \
		-o $OUTDIRARRIBA$FILE.fusions.tsv -O $OUTDIRARRIBA$FILE.fusions.discarded.tsv \
		-a $GENOMEDIR -g $ANNOTATION -b $BLACKLIST \
		-T -P > $LOGSDIR$FILE.arriba.log 2>$LOGSDIR$FILE.arriba.err
		if [ `echo $?` != 0 ]; then
			echo "An error occured during STAR/Arriba run, check the log files in $REALOUT/$FILE/LOGS/ for more details"
			exit 1
		else
			:
		fi
		wait
		if [ `echo $?` != 0 ]; then
			echo "An error occured during STAR/Arriba run, check the log files in $REALOUT/$FILE/LOGS/ for more details"
			exit 1
		else
			:
		fi
		mv $OUTDIRALIGN$FILE.Log.std.out $LOGSDIR$FILE.STAR.log
		samtools index ${OUTDIRALIGN}${FILE}.Aligned.sortedByCoord.out.bam

		if [ "$CUSTOMLIST" != "false" ]; then
			echo " Parsing custom HLA list:" `date +"%T"` | sed "s/^/[NeoFuse] /"
			python3 /usr/local/bin/source/custom_hla_parser.py --custom_list $CUSTOMLIST --output_dir ${OUTDIROPTI}"/" --sample_name $FILE
		else
			# YARA + OptiType
			## YARA
			echo " YARA Run started at:" `date +"%T"` | sed "s/^/[NeoFuse] /"
			yara_mapper --version-check 0 -e 3 -t $CORES -f bam $YARAIDX $READ1  $READ2 | \
    	        samtools view -@ $CORES -h -F 4 -b1 -o $TEMPDIROPTI"/"$FILE"_mapped_1.bam"
			# mkfifo R1 R2
			# yara_mapper --version-check 0 -e 3 -t $CORES -f bam $YARAIDX $READ1  $READ2 | \
			# samtools view -@ $CORES -h -F 4 -b1 | tee R1 R2 > /dev/null &
			# samtools view -@ $CORES -h -f 0x40 -b1 R1 > $TEMPDIROPTI"/"$FILE"_mapped_1.bam" &
			# samtools view -@ $CORES -h -f 0x80 -b1 R2 > $TEMPDIROPTI"/"$FILE"_mapped_2.bam" &
			# wait
			# rm -f R1 R2
			## Optitype
			echo " OptiType Run started at:" `date +"%T"` | sed "s/^/[NeoFuse] /"
			python $Optitype -i $TEMPDIROPTI"/"$FILE"_mapped_1.bam" \
				--rna -v -o $TEMPDIROPTI > $LOGSDIR$FILE.optitype.log 2>&1
			# python $Optitype -i $TEMPDIROPTI$FILE"_mapped_1.bam" $TEMPDIROPTI$FILE"_mapped_2.bam" -e 1 -b 0.009 -v --rna -o $TEMPDIROPTI > $LOGSDIR$FILE.optitype.log 2>&1
			if [ `echo $?` != 0 ]; then
				echo "An error occured during OptiType run, check $REALOUT/$FILE/LOGS/$FILE.optitype.log for more details"
				exit 1
			else
				:
			fi
			## Save output and remove temporary files
			OutFile1=${OUTDIROPTI}$FILE"_HLA_Optitype.txt"
			OutFile2=${OUTDIROPTI}$FILE"_coverage_plot.pdf"
			tmpFile=`ls ${TEMPDIROPTI}/*/*result.tsv`
			pdfFile=`ls ${TEMPDIROPTI}/*/*_coverage_plot.pdf`
			mv $pdfFile $OutFile2
			tail -1 $tmpFile | cut -f 2-7 | tr "\t" "\n" | sort | uniq > $OutFile1
		fi

		# Asign Reads to features (featureCounts)
		echo " featureCounts Run started at: "`date +"%T"` | sed "s/^/[NeoFuse] /"
		featureCounts -p -t exon -T $CORES \
			-g gene_name \
			-a $ANNOTATION \
			-o $OUTDIRCOUNTS$FILE.counts.txt \
			$OUTDIRALIGN$FILE.Aligned.sortedByCoord.out.bam > $LOGSDIR$FILE.featureCounts.log 2>&1

		if [ `echo $?` != 0 ]; then
			echo "An error occured during featureCounts run, check $REALOUT/$FILE/LOGS/$FILE.featureCounts.log for more details"
			exit 1
		else
			:
		fi
	else
		echo "[-------------------------------- [NeoFuse] --------------------------------]"
		echo
		echo " Single End (SE) Reads detected: commencing processing" | sed "s/^/[NeoFuse] /"
		echo " Processing file $NAME1" | sed "s/^/[NeoFuse] /"
		echo " STAR Run started at:" `date +"%T"` | sed "s/^/[NeoFuse] /"
		STAR --runThreadN $STARTHREADS \
		--runMode alignReads \
		--genomeDir $STARINDEX \
		--readFilesIn $READ1 $ZIPPED \
		--outFileNamePrefix $OUTDIRALIGN$FILE"." \
		--genomeLoad NoSharedMemory \
		--outReadsUnmapped Fastx \
		--outStd BAM_Unsorted \
		--outBAMcompression 0 \
		--outSAMtype BAM Unsorted 2>$LOGSDIR$FILE.STAR.err |
			samtools sort -@ $SAMTOOLSTHREADS -O BAM -o ${OUTDIRALIGN}${FILE}.Aligned.sortedByCoord.out.bam /dev/stdin \
			> $LOGSDIR$FILE.samtools.log 2>$LOGSDIR$FILE.samtools.err &

		# Arriba
		echo " Arriba Run started at:" `date +"%T"` | sed "s/^/[NeoFuse] /"
		STAR --runThreadN $ARRIBATHREADS \
		--genomeDir $STARINDEX --genomeLoad NoSharedMemory --limitBAMsortRAM $RAMLIMIT \
		--readFilesIn $READ1 $ZIPPED \
		--outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 \
		--outFilterMultimapNmax 1 --outFilterMismatchNmax 3 \
		--outFileNamePrefix $OUTDIRARRIBA$FILE"." \
		--chimSegmentMin 10 --chimOutType WithinBAM SoftClip --chimJunctionOverhangMin 10 --chimScoreMin 1 \
		--chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --alignSJstitchMismatchNmax 5 -1 5 5 \
		--chimSegmentReadGapMax 3 |
		arriba \
		-x /dev/stdin \
		-o $OUTDIRARRIBA$FILE.fusions.tsv -O $OUTDIRARRIBA$FILE.fusions.discarded.tsv \
		-a $GENOMEDIR -g $ANNOTATION -b $BLACKLIST \
		-T -P > $LOGSDIR$FILE.arriba.log 2>$LOGSDIR$FILE.arriba.err
		if [ `echo $?` != 0 ]; then
			echo "An error occured during STAR/Arriba run, check the log files in $REALOUT/$FILE/LOGS/ for more details"
			exit 1
		else
			:
		fi
		wait
		if [ `echo $?` != 0 ]; then
			echo "An error occured during STAR/Arriba run, check the log files in $REALOUT/$FILE/LOGS/ for more details"
			exit 1
		else
			:
		fi
		mv $OUTDIRALIGN$FILE.Log.std.out $LOGSDIR$FILE.STAR.log
		samtools index ${OUTDIRALIGN}${FILE}.Aligned.sortedByCoord.out.bam

		if [ "$CUSTOMLIST" != "false" ]; then
			echo " Parsing custom HLA list:" `date +"%T"` | sed "s/^/[NeoFuse] /"
			python3 /usr/local/bin/source/custom_hla_parser.py --custom_list $CUSTOMLIST --output_dir ${OUTDIROPTI}"/" --sample_name $FILE
		else
			# YARA + OptiType
			## YARA
			echo " YARA Run started at:" `date +"%T"` | sed "s/^/[NeoFuse] /"
			yara_mapper --version-check 0 -e 3 -t $CORES -f bam $YARAIDX $READ1 | \
    	        samtools view -@ $CORES -h -F 4 -b1 -o $TEMPDIROPTI"/"$FILE"_mapped_1.bam"
			## Optitype
			echo " OptiType Run started at:" `date +"%T"` | sed "s/^/[NeoFuse] /"
			python $Optitype -i $TEMPDIROPTI"/"$FILE"_mapped_1.bam" \
				--rna -v -o $TEMPDIROPTI > $LOGSDIR$FILE.optitype.log 2>&1
			# python $Optitype -i $TEMPDIROPTI$FILE"_mapped_1.bam" $TEMPDIROPTI$FILE"_mapped_2.bam" -e 1 -b 0.009 -v --rna -o $TEMPDIROPTI > $LOGSDIR$FILE.optitype.log 2>&1
			if [ `echo $?` != 0 ]; then
				echo "An error occured during OptiType run, check $REALOUT/$FILE/LOGS/$FILE.optitype.log for more details"
				exit 1
			else
				:
			fi
			## Save output and remove temporary files
			OutFile1=${OUTDIROPTI}$FILE"_HLA_Optitype.txt"
			OutFile2=${OUTDIROPTI}$FILE"_coverage_plot.pdf"
			tmpFile=`ls ${TEMPDIROPTI}/*/*result.tsv`
			pdfFile=`ls ${TEMPDIROPTI}/*/*_coverage_plot.pdf`
			mv $pdfFile $OutFile2
			tail -1 $tmpFile | cut -f 2-7 | tr "\t" "\n" | sort | uniq > $OutFile1
		fi
<<<<<<< HEAD
		
=======

>>>>>>> c0625b6e43efb01d4dfafc8ad3c31e04ab669e8c
		# Asign Reads to features (featureCounts)
		echo " featureCounts Run started at: "`date +"%T"` | sed "s/^/[NeoFuse] /"
		featureCounts -t exon -T $CORES \
			-g gene_name \
			-a $ANNOTATION \
			-o $OUTDIRCOUNTS$FILE.counts.txt \
			$OUTDIRALIGN$FILE.Aligned.sortedByCoord.out.bam > $LOGSDIR$FILE.featureCounts.log 2>&1
		if [ `echo $?` != 0 ]; then
			echo "An error occured during featureCounts run, check $REALOUT/$FILE/LOGS/$FILE.featureCounts.log for more details"
			exit 1
		else
			:
		fi
	fi

	# Convert raw counts to TPM/RPKM (Rscript)
	echo " Converting Raw Counts to TPM and FPKM:" `date +"%T"` | sed "s/^/[NeoFuse] /"
	Rscript /usr/local/bin/source/counts_to_tpm.R \
		$OUTDIRCOUNTS$FILE.counts.txt \
		$OUTDIRTPM$FILE.tpm.txt \
		$OUTDIRRPKM$FILE.rpkm.txt > $LOGSDIR$FILE.counts_to_tpm.log 2>&1

	if [ `echo $?` != 0 ]; then
		echo "An error occured during conversion of Raw Counts to TPM and FPKM, check $REALOUT/$FILE/LOGS/$FILE.counts_to_tpm.log for more details"
		exit 1
	else
		:
	fi
	# Cleave peptides
	echo " Searching for peptides of length $PEPLEN:" `date +"%T"` | sed "s/^/[NeoFuse] /"
	python3 /usr/local/bin/source/cleave_peptides.py -i $OUTDIRARRIBA$FILE.fusions.tsv -o $OUTDIRCLEAVEPEP$FILE -p $PEPLEN > $LOGSDIR$FILE.cleave.log 2>&1
	if [ `echo $?` != 0 ]; then
		echo "An error occured while searching for peptides, check $REALOUT/$FILE/LOGS/$FILE.cleave_peptides.log for more details"
		exit 1
	else
		:
	fi

	if [ "$NETMHCPAN" == "false" ]; then
		# MHCFlurry and final output
		## create the associations and run files
		echo " MHCFlurry Run started at: "`date +"%T"` | sed "s/^/[NeoFuse] /"
		for filename in $OUTDIRCLEAVEPEP*.tsv; do
			fou=$(echo ${filename##*/} | sed s/_xeno// | sed s/\.tsv//)
   			python3 /usr/local/bin/source/association.py -x $filename -l ${OUTDIROPTI}"/"$FILE"_HLA_Optitype.txt" -o $FINALTMP$fou -c $CORES > $LOGSDIR$FILE.association.log 2>&1
		done
		if [ `echo $?` != 0 ]; then
			echo "An error occured while creating the assosciation files, check $REALOUT/$FILE/LOGS/$FILE.association.log for more details"
			exit 1
		else
			:
		fi

		## start MHCFlurry predictions
		for i in $FINALTMP*_TEST_OUT.sh; do
			name=$(echo ${i##*/} | sed s/_TEST_OUT.sh//)
			sh $i > $LOGSDIR$name"_MHCFlurry.log" 2>&1
			if [ `echo $?` != 0 ]; then
				echo "An error occured during the MHCFlurry run, check $REALOUT/$FILE/LOGS/$name\"_MHCFlurry.log\" for more details"
				exit 1
			else
				:
			fi
		done
		wait
		sleep 30 # Extra time to release resources
		echo " Creating Final Ouptut:" `date +"%T"` | sed "s/^/[NeoFuse] /"
		for j in $PEPLEN; do
			python3 /usr/local/bin/source/build_temp.py -a $FINALTMP*$j"_ASSOCIATIONS_OUT.txt" -o $FINALTMP$FILE"_"$j > $LOGSDIR$FILE.final.log 2>&1
			if [ `echo $?` != 0 ]; then
				echo "An error occured while creating the final output files, check $REALOUT/$FILE/LOGS/$FILE.final.log for more details"
				exit 1
			else
				:
			fi
		done

		for file in $FINALTMP$FILE"_"*_unsupported.txt; do
			cat $file > $FINALOUTDIR$FILE"_unsupported.txt"
		done
	else
		# netMHCpan and final output
		## export env paths
		export TMPDIR=$FINALTMP
		export NHOME=$NETMHCPAN
		export NETMHCpan=$NETMHCPAN/Linux_x86_64
		## create the association file and run files
		echo " netMHCpan Run started at: "`date +"%T"` | sed "s/^/[NeoFuse] /"
		for filename in $OUTDIRCLEAVEPEP*.tsv; do
			fou=$(echo ${filename##*/} | sed s/_xeno// | sed s/\.tsv//)
   			python3 /usr/local/bin/source/association_netMHCpan.py -x $filename -l ${OUTDIROPTI}"/"$FILE"_HLA_Optitype.txt" -o $FINALTMP$fou -p $NETMHCPAN -c $CORES > $LOGSDIR$FILE.association.log 2>&1
		done
		if [ `echo $?` != 0 ]; then
			echo "An error occured while creating the assosciation files, check $REALOUT/$FILE/LOGS/$FILE.association.log for more details"
			exit 1
		else
			:
		fi

		## start netMHCpan predictions
		for i in $FINALTMP*_TEST_OUT.sh; do
			name=$(echo ${i##*/} | sed s/_TEST_OUT.sh//)
 			sh $i > $LOGSDIR$name"_netMHCpan.log" 2>&1
 			if [ `echo $?` != 0 ]; then
 				echo "An error occured during the netMHCpan run, check $REALOUT/$FILE/LOGS/$name\"_netMHCpan.log\" for more details"
 				exit 1
	 		else
	 			:
	 		fi
		done
		wait
		sleep 30 # Extra time to release resources
		echo " Creating Final Ouptut:" `date +"%T"` | sed "s/^/[NeoFuse] /"
		for j in $PEPLEN; do
			python3 /usr/local/bin/source/build_temp_netMHCpan.py -a $FINALTMP*$j"_ASSOCIATIONS_OUT.txt" -o $FINALTMP$FILE"_"$j > $LOGSDIR$FILE.final.log 2>&1
			if [ `echo $?` != 0 ]; then
				echo "An error occured while creating the final output files, check $REALOUT/$FILE/LOGS/$FILE.final.log for more details"
				exit 1
			else
				:
			fi
		done
	fi
<<<<<<< HEAD
	
=======

>>>>>>> c0625b6e43efb01d4dfafc8ad3c31e04ab669e8c
	echo "Fusion	Gene1	Gene2	Breakpoint1	Breakpoint2	HLA_type	Fusion_Peptide	IC50	Rank	Event_Type	Stop_Codon	Confidence" > $FINALOUTDIR$FILE"_unfiltered.tsv"
	for file in $FINALTMP$FILE*_final.tsv; do
		cat $file | sed 1d >> $FINALOUTDIR$FILE"_unfiltered.tsv"
	done

	# Filter the results file
	python3 /usr/local/bin/source/filter.py -i $FINALOUTDIR$FILE"_unfiltered.tsv" -o $FINALOUTDIR$FILE"_tmp_filtered.tsv" -t $THRESHOLD -T $RANK -c $CONFIDENCE >> $LOGSDIR$FILE.final.log 2>&1
	if [ `echo $?` != 0 ]; then
		echo "An error occured while filtering the final output files, check $REALOUT/$FILE/LOGS/$FILE.final.log for more details"
		exit 1
	else
		:
	fi

	# Add TPM and compute average TPM on final results
	python3 /usr/local/bin/source/final_out.py -t $OUTDIRTPM$FILE.tpm.txt -i $FINALOUTDIR$FILE"_tmp_filtered.tsv" -o $FINALOUTDIR$FILE"_filtered.tsv" >> $LOGSDIR$FILE.final.log 2>&1
	if [ `echo $?` != 0 ]; then
		echo "An error occured while appending TPM values to the final output files, check $REALOUT/$FILE/LOGS/$FILE.final.log for more details"
		exit 1
	else
		:
	fi
	## Clean up
	echo " Removing Intermediate Files:" `date +"%T"` | sed "s/^/[NeoFuse] /"
	if [ "$CUSTOMLIST" != "false" ]; then
		mkdir -p $OUTDIR"Custom_HLAs/"
		mv ${OUTDIROPTI}"/"$FILE"_HLA_Optitype.txt" $OUTDIR"Custom_HLAs/"$FILE"_custom_HLA_I.txt"
		mv ${OUTDIROPTI}"/"$FILE"_custom_HLA_II.txt" $OUTDIR"Custom_HLAs/"
		rm -rf ${OUTDIROPTI}
	else
		rm -rf $TEMPDIROPTI
	fi
	rm $FINALOUTDIR$FILE"_tmp_filtered.tsv"
	rm $OUTDIRARRIBA$FILE".Log.progress.out"
	rm $OUTDIRARRIBA$FILE".SJ.out.tab"
	mv $OUTDIRARRIBA$FILE".Log."* $LOGSDIR
	if [ "$KEEPBAM" == "false" ]; then
		rm -rf $OUTDIR"/STAR/"
	else
		:
	fi
	rm -rf $OUTDIRCOUNTS
	rm -rf $OUTDIRRPKM
	rm -rf $OUTDIRCLEAVEPEP
	rm -rf $FINALTMP
	echo " The run has succesfully finished at:" `date +"%T"` | sed "s/^/[NeoFuse] /"
	echo " All the results and log files can be found in $OUTDIR/$FILE/" | sed "s/^/[NeoFuse] /"
	echo "[-------------------------------- [NeoFuse] --------------------------------]"
	echo
done

rm $OUTDIR"/in.tsv"
rm $IN
