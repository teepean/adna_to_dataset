#!/bin/bash
command -v pileupCaller >/dev/null 2>&1 || { echo >&2 "pileupCaller is required.  Aborting."; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo >&2 "samtools is required.  Aborting."; exit 1; }

#!/bin/bash

# Check if hs37d5.fa exists in the reference directory
if [ ! -f reference/hs37d5.fa ]; then
    echo
    echo "You need to have hs37d5.fa in the reference directory."
    echo "Do you want to download it? (Y/N)"
    read -r choice
    case $choice in
        [Yy]* )
            echo "Downloading..."
            cd reference || exit
            curl --insecure -O https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
            echo "Decompressing..."
            bgzip -d -@ 2 hs37d5.fa.gz
            echo "Indexing hs37d5.fa..."
            samtools faidx hs37d5.fa
            cd ..
            ;;
        [Nn]* )
            echo "Skipping download of hs37d5.fa."
            ;;
        * )
            echo "Invalid response. Exiting."
            exit 1
            ;;
    esac
    echo
fi

# Check if hg19.fa exists in the reference directory
if [ ! -f reference/hg19.fa ]; then
    echo
    echo "You need to have hg19.fa in the reference directory."
    echo "Do you want to download it? (Y/N)"
    read -r choice
    case $choice in
        [Yy]* )
            echo "Downloading..."
            cd reference || exit
            curl --insecure -O https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz
            echo "Decompressing..."
            bgzip -d -@ 2 hg19.fa
            echo "Indexing hg19.fa..."
            samtools faidx hg19.fa
            cd ..
            ;;
        [Nn]* )
            echo "Skipping download of hg19.fa."
            ;;
        * )
            echo "Invalid response. Exiting."
            exit 1
            ;;
    esac
    echo
fi

# Prompt user for Population name
read -p "Enter Population name: " POPNAME

# Prompt user for Output name
read -p "Enter Output name: " OUTPUTNAME

# Navigate to the source directory and list all .bam files, storing them in bamlist1
cd source
bamlist1=($(ls *.bam *.cram))
IFS=" " bamlist_space="${bamlist1[*]}"

# Initialize array2 as an empty array
array2=()

# Loop through bamlist1, remove the .bam extension from each element, and add it to array2
for i in "${bamlist_space[@]}"
do
  array2+=("${i%.bam}")
done

# Convert array2 into a comma-separated string
array2_string=$(echo "$array2" | tr ' ' ',')

read -p "Enter your choice (hs37d5 or hg19): " choice

if [[ "$choice" == "hs37d5" ]]; then
    echo "Running command"
    samtools mpileup -B -q 30 -Q 30 -l ../positions/v42.4.1240K.pos -f ../reference/hs37d5.fa ${bamlist_space} | pileupCaller --randomHaploid --sampleNames $array2_string --samplePopName $POPNAME -f ../positions/v42.4.1240K.snp -p ../target/$OUTPUTNAME
elif [[ "$choice" == "hg19" ]]; then
    echo "Running command"
    samtools mpileup -B -q 30 -Q 30 -l ../positions/v42.4.hg19.pos -f ../reference/hg19.fa ${bamlist_space} | sed 's/chr//' | pileupCaller --randomHaploid --sampleNames $array2_string --samplePopName $POPNAME -f ../positions/v42.4.1240K.snp -p ../target/$OUTPUTNAME
else
    echo "Invalid choice"
fi

cd ..

