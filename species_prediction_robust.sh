#!/bin/bash

# ===================================================
# UNIVERSAL Automated Species Prediction Pipeline
# Works on: macOS (M1/M2/Intel) & Linux (Ubuntu/Debian/CentOS)
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

# Check if input files exist
if [ ! -f "$READ1" ] || [ ! -f "$READ2" ]; then
    echo "Error: Input files $READ1 or $READ2 not found."
    exit 1
fi

echo "==========================================="
echo "Step 1: Creating Subset (100k reads)"
echo "==========================================="
head -400000 "$READ1" > "$SUB1"
head -400000 "$READ2" > "$SUB2"

echo "==========================================="
echo "Step 2: Running SPAdes"
echo "==========================================="
# Auto-detect cores for Mac/Linux to speed up assembly
if [[ "$OSTYPE" == "darwin"* ]]; then
    CORES=$(sysctl -n hw.ncpu)
else
    CORES=$(nproc)
fi

spades.py --isolate -t "$CORES" -1 "$SUB1" -2 "$SUB2" -o "$OUTPUT_DIR"

if [ ! -f "$OUTPUT_DIR/contigs.fasta" ]; then
    echo "Assembly failed. Check SPAdes logs."
    exit 1
fi

echo "==========================================="
echo "Step 3: Extract Longest Contig"
echo "==========================================="
# This AWK block is POSIX compliant (Mac/Linux safe)
awk '
    /^>/ {if (seq) print length(seq), header, seq; header=$0; seq=""; next}
    {seq=seq""$0}
    END {print length(seq), header, seq}
' "$OUTPUT_DIR/contigs.fasta" | sort -nr | head -n 1 > temp.txt

awk '{print $2; print $3}' temp.txt > "$LONGEST"
rm temp.txt

echo "Longest contig extracted to $LONGEST"

echo "==========================================="
echo "Step 4: Submitting to NCBI BLAST"
echo "==========================================="

# Use curl's built-in urlencoding for the sequence
# Extract RID using a sed pattern that works on both BSD and GNU
RESPONSE=$(curl -s -X POST "https://blast.ncbi.nlm.nih.gov/Blast.cgi" \
    --data-urlencode "CMD=Put" \
    --data-urlencode "PROGRAM=blastn" \
    --data-urlencode "DATABASE=nt" \
    --data-urlencode "QUERY=$(cat $LONGEST)")

RID=$(echo "$RESPONSE" | grep "RID =" | head -n 1 | sed 's/.*RID = //;s/ .*//')

if [ -z "$RID" ]; then
    echo "Error: Failed to get RID from NCBI."
    exit 1
fi

echo "RID: $RID"
echo "Waiting for BLAST result (Polling every 20s)..."

COUNT=0
MAX_WAIT=900 # Increased to 15 mins for massive nt database
while true; do
    # Capture status and use tr to strip all whitespace for reliable comparison
    STATUS_RAW=$(curl -s "https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=$RID&FORMAT_OBJECT=SearchInfo")
    STATUS=$(echo "$STATUS_RAW" | grep "Status=" | cut -d= -f2 | tr -d '[:space:]')

    if [[ "$STATUS" == "READY" ]]; then
        echo "Results are ready!"
        break
    elif [[ "$STATUS" == "FAILED" ]]; then
        echo "BLAST search failed at NCBI."
        exit 1
    fi

    echo "...status: $STATUS (${COUNT}s elapsed)"
    sleep 20
    COUNT=$((COUNT+20))

    if [ $COUNT -gt $MAX_WAIT ]; then
        echo "Timeout: Search took longer than 15 minutes."
        exit 1
    fi
done

echo "==========================================="
echo "Step 5: Extracting Species"
echo "==========================================="

curl -s "https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=$RID&FORMAT_TYPE=Text" > blast_result.txt

# Extract the species name from the top alignment line
# Works by looking for the first non-comment line after the alignment header
SPECIES=$(grep -A 10 "Sequences producing significant alignments" blast_result.txt \
          | grep -v "Sequences producing" \
          | grep -v "\-\-\-\-\-" \
          | grep -v "^$" \
          | head -n 1 \
          | sed 's/  */ /g' \
          | cut -d ' ' -f 1-8)

echo "-------------------------------------------"
echo "Final Predicted Species (Top Hit):"
echo "$SPECIES"
echo "-------------------------------------------"