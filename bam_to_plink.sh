#!/bin/bash

# Enable case-insensitive pattern matching
shopt -s nocasematch

# Check if hs37d5.fa exists in the reference directory
if [ ! -f reference/hs37d5.fa ]; then
  echo "You need to have hs37d5.fa in the reference directory."
  read -p "Do you want to download it? (y/n) " yn
  case $yn in
    [Yy]* )
      echo "Downloading..."
      cd reference
      curl --insecure -O https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
      echo "Decompressing..."
      bgzip -d -@ 2 hs37d5.fa.gz
      echo "Indexing hs37d5.fa..."
      samtools faidx hs37d5.fa
      cd ..
      ;;
    * ) ;;
  esac
fi

# Check if hg19.fa exists in the reference directory
if [ ! -f reference/hg19.fa ]; then
  echo "You need to have hg19.fa in the reference directory."
  read -p "Do you want to download it? (y/n) " yn
  case $yn in
    [Yy]* )
      echo "Downloading..."
      cd reference
      curl --insecure -O https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz
      echo "Decompressing..."
      bgzip -d -@ 2 hg19.fa.gz
      echo "Indexing hg19.fa..."
      samtools faidx hg19.fa
      cd ..
      ;;
    * ) ;;
  esac
fi

# Prompt user for Population name
read -p "Enter Population name: " POPNAME

# Prompt user for Output name
read -p "Enter Output name: " OUTPUTNAME

# Navigate to the source directory and list all .bam and .cram files
cd source
# Enable nullglob to ensure unmatched patterns are ignored
shopt -s nullglob

bamlist1=""
headerCheck="false"

for file in *.bam *.cram; do
  if [ "$headerCheck" = "false" ]; then
    headerCheck="true"
    HEADER_TEMP="header_temp.txt"
    samtools view -H "$file" > "$HEADER_TEMP"
    FLAG="NEITHER"

    if grep -q -b -f ../patternhs37d5.txt "$HEADER_TEMP"; then
      FLAG="HS37D5"
    fi
    if grep -q -b -f ../patternhg19.txt "$HEADER_TEMP"; then
      FLAG="HG19"
    fi
    rm "$HEADER_TEMP"
  fi
  bamlist1="$bamlist1 $file"
done

# Disable nullglob if you don't want it to affect subsequent scripts or commands
shopt -u nullglob

array2=()

for i in $bamlist1; do
  name=$(basename "$i" .bam)
  array2+=("$name")
done

array2_string=$(IFS=,; echo "${array2[*]}")

if [ "$FLAG" == "HS37D5" ]; then
  echo "Running command"
  samtools mpileup -B -q 30 -Q 30 -l ../positions/v42.4.1240K.pos -f ../reference/hs37d5.fa $bamlist1 | pileupCaller --randomHaploid --sampleNames "$array2_string" --samplePopName "$POPNAME" -f ../positions/v42.4.1240K.snp -p ../target/"$OUTPUTNAME"
elif [ "$FLAG" == "HG19" ]; then
  echo "Running command"
  samtools mpileup -B -q 30 -Q 30 -l ../positions/v42.4.hg19.pos -f ../reference/hg19.fa $bamlist1 | sed "s/chr//" | pileupCaller --randomHaploid --sampleNames "$array2_string" --samplePopName "$POPNAME" -f ../positions/v42.4.1240K.snp -p ../target/"$OUTPUTNAME"
else
  echo "Reference is not compatible"
  read -p "Press enter to continue..."
fi
