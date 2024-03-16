#!/bin/bash
command -v pileupCaller >/dev/null 2>&1 || { echo >&2 "pileupCaller is required.  Aborting."; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo >&2 "samtools is required.  Aborting."; exit 1; }

# Prompt user for Population name
read -p "Enter Population name: " POPNAME

# Prompt user for Output name
read -p "Enter Output name: " OUTPUTNAME

# Navigate to the source directory and list all .bam files, storing them in bamlist1
cd source
bamlist1=($(ls *.bam))
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
    samtools mpileup -B -q 30 -Q 30 -l ../positions/v42.4.1240K.pos -f ../reference/hs37d5.fa.gz ${bamlist_space} | pileupCaller --randomHaploid --sampleNames $array2_string --samplePopName $POPNAME -f ../positions/v42.4.1240K.snp -p ../target/$OUTPUTNAME
elif [[ "$choice" == "hg19" ]]; then
    echo "Running command"
    samtools mpileup -B -q 30 -Q 30 -l ../positions/v42.4.hg19.pos -f ../reference/hg19.fa.gz ${bamlist_space} | sed 's/chr//' | pileupCaller --randomHaploid --sampleNames $array2_string --samplePopName $POPNAME -f ../positions/v42.4.1240K.snp -p ../target/$OUTPUTNAME
else
    echo "Invalid choice"
fi

# Use the variables in the samtools command
# samtools mpileup -B -q 30 -Q 30 -l positions/v42.4.1240K.pos -f reference/hs37d5.fa source/${bamlist_space[*]} | pileupCaller --randomHaploid --sampleNames $array2_string --samplePopName $POPNAME -f positions/v42.4.1240K.snp -p target/$OUTPUTNAME

cd ..

