#!/bin/bash

# ===================================================
# UNIVERSAL Automated Species Prediction Pipeline (Single-End)
# Works on: macOS & Linux
# ===================================================

if [ -z "$1" ]; then
    echo "Usage: $0 <reads.fastq>"
    exit 1
fi

READ=$1
SUB="subset.fastq"
OUTPUT_DIR="spades_output_single"
LONGEST="longest_contig.fasta"
SHORT_QUERY="short_query.fasta"

# Check if input file exists
if [ ! -f "$READ" ]; then
    echo "Error: Input file $READ not found."
    exit 1
fi

echo "==========================================="
echo "Step 1: Creating Subset (100k reads)"
echo "==========================================="
head -400000 "$READ" > "$SUB"

echo "==========================================="
echo "Step 2: Running SPAdes (Single-End)"
echo "==========================================="

# Auto-detect cores for Mac/Linux
if [[ "$OSTYPE" == "darwin"* ]]; then
    CORES=$(sysctl -n hw.ncpu)
else
    CORES=$(nproc)
fi

spades.py --isolate -t "$CORES" -s "$SUB" -o "$OUTPUT_DIR"

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
' "$OUTPUT_DIR/contigs.fasta" | sort -nr | head -n 1 > temp.txt

awk '{print $2; print $3}' temp.txt > "$LONGEST"
rm temp.txt

echo "Longest contig extracted."

echo "==========================================="
echo "Step 4: Preparing Short BLAST Query (1000 bp)"
echo "==========================================="

HEADER=$(head -n 1 "$LONGEST")
# Extract first 1000bp safely across platforms
SEQUENCE=$(tail -n +2 "$LONGEST" | tr -d '\n' | cut -c1-1000)

echo "$HEADER" > "$SHORT_QUERY"
echo "$SEQUENCE" >> "$SHORT_QUERY"

echo "Submitting BLAST query..."

# Use curl's data-urlencode to avoid sed/newline encoding issues
RESPONSE=$(curl -s -X POST "https://blast.ncbi.nlm.nih.gov/Blast.cgi" \
    --data-urlencode "CMD=Put" \
    --data-urlencode "PROGRAM=blastn" \
    --data-urlencode "DATABASE=nt" \
    --data-urlencode "ENTREZ_QUERY=Mammalia[Organism]" \
    --data-urlencode "QUERY=$(cat $SHORT_QUERY)")

# Mac-safe RID extraction
RID=$(echo "$RESPONSE" | grep "RID =" | head -n 1 | sed 's/.*RID = //;s/ .*//')

if [ -z "$RID" ]; then
    echo "Failed to get RID from NCBI."
    echo "Server Response Snippet:"
    echo "$RESPONSE" | head -n 15
    exit 1
fi

echo "RID: $RID"
echo "Waiting for BLAST result (Polling every 20s)..."

COUNT=0
MAX_WAIT=600
while true; do
    # Capture status and strip whitespace/trailing \r for reliable comparison
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
        echo "Timeout exceeded (10 mins)."
        exit 1
    fi
done

echo "Fetching result..."
curl -s "https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=$RID&FORMAT_TYPE=Text" > blast_result.txt

echo "==========================================="
echo "Step 5: Extracting Species"
echo "==========================================="

# More robust species extraction for Mac/Linux grep versions
SPECIES=$(grep -A 10 "Sequences producing significant alignments" blast_result.txt \
          | grep -v "Sequences producing" \
          | grep -v "\-\-\-\-\-" \
          | grep -v "^$" \
          | head -n 1 \
          | sed 's/  */ /g' \
          | cut -d ' ' -f 1-8)

echo "==========================================="
echo "Final Predicted Species:"
echo "$SPECIES"
echo "==========================================="