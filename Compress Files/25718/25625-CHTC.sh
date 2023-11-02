#! /bin/bash

usage () { echo "Usage : $0 \
-s <sample name> \
-p <path to reads on transfer server> \
-r <R1 FASTQ file> \
-t <R2 FASTQ file> \
-f <Reference FASTA file for mapping> \
-e <experiment number>
"; }

while getopts s:p:r:t:f:e: opt ; do
   case $opt in
      s) SAMPLE_NAME=$OPTARG ;;
      p) TRANSFER_FASTQ_PATH=$OPTARG ;;
      r) R1_FASTQ=$OPTARG ;;
      t) R2_FASTQ=$OPTARG ;;
      f) REF_FASTA=$OPTARG ;;
      e) EXPERIMENT_NUMBER=$OPTARG ;;
      *) usage; exit 1;;
   esac
done

# error if missing parameter and exit

if [ ! "$SAMPLE_NAME" ] || [ ! "$TRANSFER_FASTQ_PATH" ] || [ ! "$R1_FASTQ" ] || [ ! "$R2_FASTQ" ] || [ ! "$REF_FASTA" ] || [ ! "$EXPERIMENT_NUMBER" ] 
then
    usage
    exit
fi

# copy FASTQ files to execute node from transfer server
cp $TRANSFER_FASTQ_PATH/$R1_FASTQ .
cp $TRANSFER_FASTQ_PATH/$R2_FASTQ .

# invoke snakemake
snakemake --cores 1 --snakefile 25625-CHTC.snakefile --config experiment_number="$EXPERIMENT_NUMBER" sample_name="$SAMPLE_NAME" r1_fastq="$R1_FASTQ" r2_fastq="$R2_FASTQ" ref_fasta="$REF_FASTA"

rm $R1_FASTQ
rm $R2_FASTQ