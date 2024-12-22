# Enable strict mode
Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

# Function to check command length
function Test-CommandLength {
   param (
       [string]$Command,
       [int]$Limit = 8000  # Setting slightly below 8,191 for safety margin
   )
   if ($Command.Length -gt $Limit) {
       Write-Warning "Command length ($($Command.Length) characters) exceeds safe limit of $Limit characters."
       Write-Warning "This may cause the command to fail."
       Write-Warning "Consider processing fewer files at once."
       
       # Optional: Ask user whether to proceed
       $proceed = Read-Host "Do you want to proceed anyway? (Y/N)"
       if ($proceed -ne 'Y') {
           return $false
       }
   }
   return $true
}

# Function to process files in batches
function Process-FilesInBatch {
   param (
       [array]$FileList,
       [string]$ReferenceType,
       [string]$PopName,
       [string]$OutputName,
       [int]$BatchSize = 10
   )

   $totalBatches = [Math]::Ceiling($FileList.Count / $BatchSize)
   $currentBatch = 1

   for ($i = 0; $i -lt $FileList.Count; $i += $BatchSize) {
       $batch = $FileList[$i..([Math]::Min($i + $BatchSize - 1, $FileList.Count - 1))]
       Write-Host "Processing batch $currentBatch of $totalBatches"

       $bamlist_space = $batch -join ' '
       $array2_string = ($batch | ForEach-Object { [System.IO.Path]::GetFileNameWithoutExtension($_) }) -join ','
       
       try {
           if ($ReferenceType -eq "HS37D5") {
               # Create a batch file to handle the pipeline
               $batchContent = @"
@echo off
..\winbin\samtools mpileup -B -q 30 -Q 30 -l ../positions/v42.4.1240K.pos -f ../reference/hs37d5.fa $bamlist_space | ^
..\winbin\pileupCaller --majorityCall --sampleNames $array2_string --samplePopName $PopName -f ../positions/v42.4.1240K.snp -p ../target/$($OutputName)_batch$currentBatch > ../target/$($OutputName)_batch$currentBatch.stats.txt 2>&1
"@
           }
           else {
               # Create a batch file to handle the pipeline
               $batchContent = @"
@echo off
..\winbin\samtools mpileup -B -q 30 -Q 30 -l ../positions/v42.4.hg19.pos -f ../reference/hg19.fa $bamlist_space | ^
..\winbin\sed "s/chr//" | ^
..\winbin\pileupCaller --majorityCall --sampleNames $array2_string --samplePopName $PopName -f ../positions/v42.4.1240K.snp -p ../target/$($OutputName)_batch$currentBatch > ../target/$($OutputName)_batch$currentBatch.stats.txt 2>&1
"@
           }

           # Write and execute the batch file
           $tempBatchPath = [System.IO.Path]::GetTempFileName() + ".bat"
           Set-Content -Path $tempBatchPath -Value $batchContent -Encoding ASCII

           # Execute the batch file
           & $tempBatchPath
           
           # Clean up
           Remove-Item $tempBatchPath -ErrorAction SilentlyContinue
       }
       catch {
           Write-Error "Batch $currentBatch failed. Error details: $_"
           return $false
       }

       $currentBatch++
   }
   return $true
}

# Check for hs37d5.fa
if (-not (Test-Path "reference\hs37d5.fa")) {
   Write-Host "`nYou need to have hs37d5.fa in the reference directory"
   Write-Host "Do you want to download it?"
   $response = Read-Host "Y/N"
   if ($response -eq 'Y') {
       Write-Host "Downloading"
       Push-Location reference
       & ..\winbin\curl.exe --insecure -O "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
       Write-Host "Decompressing"
       & ..\winbin\bgzip -d -@2 hs37d5.fa.gz
       Write-Host "Indexing hs37d5.fa"
       & ..\winbin\samtools faidx hs37d5.fa
       Pop-Location
   }
   else {
       exit
   }
   Write-Host "`n"
}

# Check for hg19.fa
if (-not (Test-Path "reference\hg19.fa")) {
   Write-Host "`nYou need to have hg19.fa in the reference directory"
   Write-Host "Do you want to download it?"
   $response = Read-Host "Y/N"
   if ($response -eq 'Y') {
       Write-Host "Downloading"
       Push-Location reference
       & ..\winbin\curl.exe --insecure -O "https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz"
       Write-Host "Decompressing"
       & ..\winbin\bgzip -d -@2 hg19.fa.gz
       Write-Host "Indexing hg19.fa"
       & ..\winbin\samtools faidx hg19.fa
       Pop-Location
   }
   else {
       exit
   }
   Write-Host "`n"
}

# Prompt for Population and Output names
$POPNAME = Read-Host "Enter Population name"
$OUTPUTNAME = Read-Host "Enter Output name"

# Navigate to source directory and process BAM/CRAM files
Push-Location source
$headerCheck = $null
$FLAG = "NEITHER"
$bamlist1 = @()

# Fixed Get-ChildItem command to properly handle multiple file types
Get-ChildItem -Path * -Include "*.bam", "*.cram" | ForEach-Object {
   if ($null -eq $headerCheck) {
       $headerCheck = $_.Name
       $HEADER_TEMP = "header_temp.txt"
       # Extract header to temp file
       & ..\winbin\samtools view -H $headerCheck > $HEADER_TEMP
       
       # Check patterns
       if (Select-String -Path $HEADER_TEMP -Pattern (Get-Content "..\patternhs37d5.txt")) {
           $FLAG = "HS37D5"
       }
       elseif (Select-String -Path $HEADER_TEMP -Pattern (Get-Content "..\patternhg19.txt")) {
           $FLAG = "HG19"
       }
       Remove-Item $HEADER_TEMP
   }
   $bamlist1 += $_.Name
}

if ($FLAG -eq "NEITHER") {
   Write-Host "Reference is not compatible"
   Read-Host "Press Enter to continue"
   Pop-Location
   exit 1
}

# Try processing all files at once first
$bamlist_space = $bamlist1 -join ' '
$array2_string = ($bamlist1 | ForEach-Object { [System.IO.Path]::GetFileNameWithoutExtension($_) }) -join ','

Write-Host "Attempting to process all files at once..."

try {
   if ($FLAG -eq "HS37D5") {
       # Create a batch file to handle the pipeline
       $batchContent = @"
@echo off
..\winbin\samtools mpileup -B -q 30 -Q 30 -l ../positions/v42.4.1240K.pos -f ../reference/hs37d5.fa $bamlist_space | ^
..\winbin\pileupCaller --majorityCall --sampleNames $array2_string --samplePopName $POPNAME -f ../positions/v42.4.1240K.snp -p ../target/$OUTPUTNAME > ../target/$OUTPUTNAME.stats.txt 2>&1
"@
   }
   else {
       # Create a batch file to handle the pipeline
       $batchContent = @"
@echo off
..\winbin\samtools mpileup -B -q 30 -Q 30 -l ../positions/v42.4.hg19.pos -f ../reference/hg19.fa $bamlist_space | ^
..\winbin\sed "s/chr//" | ^
..\winbin\pileupCaller --majorityCall --sampleNames $array2_string --samplePopName $POPNAME -f ../positions/v42.4.1240K.snp -p ../target/$OUTPUTNAME > ../target/$OUTPUTNAME.stats.txt 2>&1
"@
   }

   # Write and execute the batch file
   $tempBatchPath = [System.IO.Path]::GetTempFileName() + ".bat"
   Set-Content -Path $tempBatchPath -Value $batchContent -Encoding ASCII

   # Execute the batch file
   & $tempBatchPath
   
   # Clean up
   Remove-Item $tempBatchPath -ErrorAction SilentlyContinue
}
catch {
   Write-Error "Command failed. Error details: $_"
   Write-Host "Attempting batch processing instead..."
   $batchSuccess = Process-FilesInBatch -FileList $bamlist1 -ReferenceType $FLAG -PopName $POPNAME -OutputName $OUTPUTNAME
   if (-not $batchSuccess) {
       Write-Error "Both single command and batch processing failed"
       Pop-Location
       exit 1
   }
}

Pop-Location
Write-Host "Processing completed successfully"
Read-Host "Press Enter to continue"