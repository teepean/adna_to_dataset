#!/bin/bash

# Store the base directory (where the script is run from)
BASEDIR="$(pwd)"

# Default values
SOURCEDIR="source"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --source)
            SOURCEDIR="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 [--source /path/to/folder]"
            exit 1
            ;;
    esac
done

# Detect sort command - prefer GNU sort over uutils sort (which has bugs with large pileup data)
if command -v gnusort &> /dev/null; then
    SORT_CMD="gnusort"
elif sort --version 2>&1 | grep -q "GNU coreutils"; then
    SORT_CMD="sort"
else
    echo ""
    echo "WARNING: GNU sort not found. Your system appears to use uutils coreutils sort,"
    echo "which has a known bug with large pileup data. Please install GNU coreutils:"
    echo "  sudo apt install gnu-coreutils"
    echo ""
    read -p "Continue anyway? (Y/N): " choice
    case "$choice" in
        [Yy]* ) SORT_CMD="sort" ;;
        * ) exit 1 ;;
    esac
fi

# Check if hs37d5.fa exists
if [ ! -f reference/hs37d5.fa ]; then
    echo ""
    echo "You need to have hs37d5.fa in the reference directory"
    read -p "Do you want to download it? (Y/N): " choice
    case "$choice" in
        [Yy]* )
            echo "Downloading"
            cd reference || exit
            curl --insecure -O https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
            echo "Decompressing"
            bgzip -d -@2 hs37d5.fa.gz
            echo "Indexing hs37d5.fa"
            samtools faidx hs37d5.fa
            cd ..
            ;;
        [Nn]* )
            exit 0
            ;;
        * )
            echo "Invalid choice. Exiting."
            exit 1
            ;;
    esac
    echo ""
fi

# Check if hg19.fa exists
if [ ! -f reference/hg19.fa ]; then
    echo ""
    echo "You need to have hg19.fa in the reference directory"
    read -p "Do you want to download it? (Y/N): " choice
    case "$choice" in
        [Yy]* )
            echo "Downloading"
            cd reference || exit
            curl --insecure -O https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz
            echo "Decompressing"
            bgzip -d -@2 hg19.fa.gz
            echo "Indexing hg19.fa"
            samtools faidx hg19.fa
            cd ..
            ;;
        [Nn]* )
            exit 0
            ;;
        * )
            echo "Invalid choice. Exiting."
            exit 1
            ;;
    esac
    echo ""
fi

# Prompt user to choose positions set
echo ""
echo "Choose positions set:"
echo "  1) v42.4.1240K (current)"
echo "  2) v66.2M.aadr (new)"
read -p "Enter choice [1-2]: " POSCHOICE
case "$POSCHOICE" in
    1) POSNAME="v42.4.1240K" ;;
    2) POSNAME="v66.2M.aadr" ;;
    *) echo "Invalid choice. Exiting."; exit 1 ;;
esac

POS_HS37D5="$BASEDIR/positions/$POSNAME.pos"
POS_HG19="$BASEDIR/positions/$POSNAME.hg19.pos"
SNP_FILE="$BASEDIR/positions/$POSNAME.snp"

# Fall back to legacy snp filename for v42.4 (kept as v42.4.1240K.snp)
if [ ! -f "$SNP_FILE" ] && [ "$POSNAME" = "v42.4.1240K" ]; then
    SNP_FILE="$BASEDIR/positions/v42.4.1240K.snp"
fi

# Prompt user for Population name
read -p "Enter Population name: " POPNAME

# Prompt user for Output name
read -p "Enter Output name: " OUTPUTNAME

# Navigate to the source directory
cd "$SOURCEDIR" || exit

# Find all .bam and .cram files
bamlist=()
headerCheck=""

for file in *.bam *.cram; do
    # Skip if the glob didn't match anything
    [ -e "$file" ] || continue

    # Check header of first file only
    if [ -z "$headerCheck" ]; then
        headerCheck="$file"
        HEADER_TEMP="header_temp.txt"

        # Extract header to temp file
        samtools view -H "$headerCheck" > "$HEADER_TEMP"

        # Initialize flag
        FLAG="NEITHER"

        # Check patterns
        if grep -q -f "$BASEDIR/patternhs37d5.txt" "$HEADER_TEMP"; then
            FLAG="HS37D5"
        elif grep -q -f "$BASEDIR/patternhg19.txt" "$HEADER_TEMP"; then
            FLAG="HG19"
        fi

        rm "$HEADER_TEMP"
    fi

    bamlist+=("$file")
done

# Create array of sample names (without extensions)
array2=()
for file in "${bamlist[@]}"; do
    name="${file%.*}"
    array2+=("$name")
done

# Join array2 with commas
array2_string=$(IFS=,; echo "${array2[*]}")

# Run appropriate command based on reference genome
if [ "$FLAG" = "HS37D5" ]; then
    # Legacy filename for v42.4 hs37d5 positions
    if [ "$POSNAME" = "v42.4.1240K" ]; then
        POS_FILE="$BASEDIR/positions/v42.4.1240K.pos"
    else
        POS_FILE="$POS_HS37D5"
    fi
    echo "Running command"
    samtools mpileup -B -q 30 -Q 20 -l "$POS_FILE" -f "$BASEDIR/reference/hs37d5.fa" "${bamlist[@]}" | \
        awk 'BEGIN{OFS="\t"} {if($1=="X")$1=23; else if($1=="Y")$1=24; else if($1=="MT")$1=90; print}' | \
        $SORT_CMD -t "$(printf '\t')" -k1,1n -k2,2n | \
        pileupCaller --randomHaploid --sampleNames "$array2_string" --samplePopName "$POPNAME" \
        -f "$SNP_FILE" -p "$BASEDIR/target/$OUTPUTNAME" > "$BASEDIR/target/$OUTPUTNAME.stats.txt" 2>&1
elif [ "$FLAG" = "HG19" ]; then
    # Legacy filename for v42.4 hg19 positions
    if [ "$POSNAME" = "v42.4.1240K" ]; then
        POS_FILE="$BASEDIR/positions/v42.4.hg19.pos"
    else
        POS_FILE="$POS_HG19"
        # Auto-generate chr-prefixed positions file from numeric one if missing
        if [ ! -f "$POS_FILE" ] && [ -f "$BASEDIR/positions/$POSNAME.pos" ]; then
            echo "Generating $POS_FILE from $POSNAME.pos"
            awk '{print "chr"$0}' "$BASEDIR/positions/$POSNAME.pos" > "$POS_FILE"
        fi
    fi
    echo "Running command"
    samtools mpileup -B -q 30 -Q 20 -l "$POS_FILE" -f "$BASEDIR/reference/hg19.fa" "${bamlist[@]}" | \
        sed "s/chr//" | \
        awk 'BEGIN{OFS="\t"} {if($1=="X")$1=23; else if($1=="Y")$1=24; else if($1=="MT")$1=90; print}' | \
        $SORT_CMD -t "$(printf '\t')" -k1,1n -k2,2n | \
        pileupCaller --randomHaploid --sampleNames "$array2_string" --samplePopName "$POPNAME" \
        -f "$SNP_FILE" -p "$BASEDIR/target/$OUTPUTNAME" > "$BASEDIR/target/$OUTPUTNAME.stats.txt" 2>&1
else
    echo "Reference is not compatible"
    exit 1
fi

cd ..

echo "Done! Press Enter to continue..."
read -r
