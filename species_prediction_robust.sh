#!/bin/bash

# ===================================================
# FINAL Automated Species Prediction Pipeline
# Subset + SPAdes + NCBI HTTP BLAST
# ===================================================

if [ -z "$1" ]; then
    echo "Usage: $0 <ERR_base_name>"
    exit 1
fi

BASE=$1
READ1="${BASE}_1.fastq"
READ2="${BASE}_2.fastq"

SUB1="subset_1.fastq"
SUB2="subset_2.fastq"

OUTPUT_DIR="spades_output"
LONGEST="longest_contig.fasta"

echo "==========================================="
echo "Step 1: Creating Subset (100k reads)"
echo "==========================================="

head -400000 $READ1 > $SUB1
head -400000 $READ2 > $SUB2

echo "==========================================="
echo "Step 2: Running SPAdes"
echo "==========================================="

spades.py --isolate -1 $SUB1 -2 $SUB2 -o $OUTPUT_DIR

if [ ! -f "$OUTPUT_DIR/contigs.fasta" ]; then
    echo "Assembly failed."
    exit 1
fi

echo "==========================================="
echo "Step 3: Extract Longest Contig"
echo "==========================================="

awk '
    /^>/ {if (seq) print length(seq), header, seq; header=$0; seq=""; next}
    {seq=seq""$0}
    END {print length(seq), header, seq}
' $OUTPUT_DIR/contigs.fasta | sort -nr | head -n 1 > temp.txt

> $LONGEST
awk '{print $2 >> "'$LONGEST'"; print $3 >> "'$LONGEST'"}' temp.txt
rm temp.txt

echo "Longest contig extracted."

echo "==========================================="
echo "Step 4: Submitting to NCBI BLAST"
echo "==========================================="

# URL encode FASTA properly
QUERY=$(cat $LONGEST | sed ':a;N;$!ba;s/\n/%0A/g')

RID=$(curl -s \
  -d "CMD=Put&PROGRAM=blastn&DATABASE=nt&QUERY=$QUERY" \
  https://blast.ncbi.nlm.nih.gov/Blast.cgi | \
  grep RID | awk '{print $3}')

if [ -z "$RID" ]; then
    echo "Failed to get RID."
    exit 1
fi

echo "RID: $RID"
echo "Waiting for BLAST result..."

# Poll with timeout (max 5 minutes)
COUNT=0
while true; do
    STATUS=$(curl -s \
      "https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=$RID&FORMAT_OBJECT=SearchInfo" | \
      grep Status= | cut -d= -f2)

    if [ "$STATUS" == "READY" ]; then
        break
    fi

    sleep 5
    COUNT=$((COUNT+5))

    if [ $COUNT -gt 300 ]; then
        echo "Timeout exceeded."
        exit 1
    fi
done

echo "Fetching result..."

curl -s \
"https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=$RID&FORMAT_TYPE=Text" \
> blast_result.txt

echo "==========================================="
echo "Step 5: Extracting Species"
echo "==========================================="

SPECIES=$(grep -A1 "Sequences producing significant alignments" blast_result.txt \
          | tail -n 1 \
          | awk '{$1=""; print $0}' )

echo "==========================================="
echo "Final Predicted Species:"
echo "$SPECIES"
echo "==========================================="
