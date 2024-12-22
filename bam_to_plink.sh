#!/bin/bash

# Check for required programs
for prog in samtools pileupCaller; do
    if ! command -v $prog &> /dev/null; then
        echo "Error: $prog is not installed or not in PATH"
        echo "Please install $prog and try again"
        exit 1
    fi
done

# Check for hs37d5.fa
if [ ! -f "reference/hs37d5.fa" ]; then
    echo
    echo "You need to have hs37d5.fa in the reference directory"
    echo "Do you want to download it? (y/n)"
    read -n 1 answer
    echo
    
    if [ "$answer" = "y" ] || [ "$answer" = "Y" ]; then
        echo "Downloading"
        cd reference
        curl --insecure -O https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
        echo "Decompressing"
        bgzip -d -@ 2 hs37d5.fa.gz
        echo "Indexing hs37d5.fa"
        samtools faidx hs37d5.fa
        cd ..
    else
        exit 0
    fi
    echo
fi

# Check for hg19.fa
if [ ! -f "reference/hg19.fa" ]; then
    echo
    echo "You need to have hg19.fa in the reference directory"
    echo "Do you want to download it? (y/n)"
    read -n 1 answer
    echo
    
    if [ "$answer" = "y" ] || [ "$answer" = "Y" ]; then
        echo "Downloading"
        cd reference
        curl --insecure -O https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz
        echo "Decompressing"
        bgzip -d -@ 2 hg19.fa.gz
        echo "Indexing hg19.fa"
        samtools faidx hg19.fa
        cd ..
    else
        exit 0
    fi
    echo
fi

# Prompt for population name and output name
read -p "Enter Population name: " POPNAME
read -p "Enter Output name: " OUTPUTNAME

# Navigate to source directory
cd source

# Initialize variables
headerCheck=""
FLAG="NEITHER"
bamlist1=""

# Process BAM/CRAM files
for file in *.bam *.cram; do
    # Skip if no files found
    [[ -e "$file" ]] || continue
    
    bamlist1="$bamlist1 $file"
    
    # Process header for the first file only
    if [ -z "$headerCheck" ]; then
        headerCheck=$file
        HEADER_TEMP="header_temp.txt"
        samtools view -H "$headerCheck" > "$HEADER_TEMP"
        
        # Check patterns
        if grep -f ../patternhs37d5.txt "$HEADER_TEMP" > /dev/null; then
            FLAG="HS37D5"
        elif grep -f ../patternhg19.txt "$HEADER_TEMP" > /dev/null; then
            FLAG="HG19"
        fi
        rm "$HEADER_TEMP"
    fi
done

# Create sample names list
array2=""
for file in $bamlist1; do
    name=$(basename "$file" | sed 's/\.[^.]*$//')
    array2="${array2},${name}"
done
array2_string=${array2:1}  # Remove first comma

# Run appropriate command based on reference type
if [ "$FLAG" = "HS37D5" ]; then
    echo "Running command"
    samtools mpileup -B -q 30 -Q 30 -l ../positions/v42.4.1240K.pos -f ../reference/hs37d5.fa $bamlist1 | \
    pileupCaller --majorityCall --sampleNames "$array2_string" --samplePopName "$POPNAME" \
    -f ../positions/v42.4.1240K.snp -p "../target/$OUTPUTNAME" > "../target/$OUTPUTNAME.stats.txt" 2>&1
elif [ "$FLAG" = "HG19" ]; then
    echo "Running command"
    samtools mpileup -B -q 30 -Q 30 -l ../positions/v42.4.hg19.pos -f ../reference/hg19.fa $bamlist1 | \
    sed "s/chr//" | \
    pileupCaller --majorityCall --sampleNames "$array2_string" --samplePopName "$POPNAME" \
    -f ../positions/v42.4.1240K.snp -p "../target/$OUTPUTNAME" > "../target/$OUTPUTNAME.stats.txt" 2>&1
else
    echo "Reference is not compatible"
    exit 1
fi
