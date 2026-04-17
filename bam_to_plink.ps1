param(
    [string]$Source = "source"
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

# Base directory (where the script was launched from)
$BASEDIR = (Get-Location).Path

# Tool paths
$SAMTOOLS     = Join-Path $BASEDIR "winbin\samtools.exe"
$BGZIP        = Join-Path $BASEDIR "winbin\bgzip.exe"
$CURL         = Join-Path $BASEDIR "winbin\curl.exe"
$SED          = Join-Path $BASEDIR "winbin\sed.exe"
$GAWK         = Join-Path $BASEDIR "winbin\gawk\gawk.exe"
$USORT        = Join-Path $BASEDIR "winbin\usort.exe"
$PILEUPCALLER = Join-Path $BASEDIR "winbin\pileupCaller.exe"

# Function to check command length for cmd.exe's 8191-char limit
function Test-CommandLength {
    param (
        [string]$Command,
        [int]$Limit = 8000
    )
    if ($Command.Length -gt $Limit) {
        Write-Warning "Command length ($($Command.Length) characters) exceeds safe limit of $Limit characters."
        Write-Warning "This may cause the command to fail."
        Write-Warning "Consider processing fewer files at once."
        $proceed = Read-Host "Do you want to proceed anyway? (Y/N)"
        if ($proceed -ne 'Y') {
            return $false
        }
    }
    return $true
}

function Invoke-Pipeline {
    param (
        [array]$FileList,
        [string]$ReferenceType,
        [string]$PopName,
        [string]$OutputName,
        [string]$BatchSuffix = ""
    )

    $bamlist_space = $FileList -join ' '
    $array2_string = ($FileList | ForEach-Object { [System.IO.Path]::GetFileNameWithoutExtension($_) }) -join ','

    $outPath   = Join-Path $BASEDIR "target\$OutputName$BatchSuffix"
    $statsPath = Join-Path $BASEDIR "target\$OutputName$BatchSuffix.stats.txt"
    $snpPath   = Join-Path $BASEDIR "positions\v42.4.1240K.snp"
    $gawkExpr  = 'BEGIN{OFS=\"\t\"} {if($1==\"X\")$1=23; else if($1==\"Y\")$1=24; else if($1==\"MT\")$1=90; print}'
    $tab       = "`t"

    if ($ReferenceType -eq "HS37D5") {
        $posPath = Join-Path $BASEDIR "positions\v42.4.1240K.pos"
        $refPath = Join-Path $BASEDIR "reference\hs37d5.fa"
        $prefix  = "`"$SAMTOOLS`" mpileup -B -q 30 -Q 20 -l `"$posPath`" -f `"$refPath`" $bamlist_space"
    } else {
        $posPath = Join-Path $BASEDIR "positions\v42.4.hg19.pos"
        $refPath = Join-Path $BASEDIR "reference\hg19.fa"
        $prefix  = "`"$SAMTOOLS`" mpileup -B -q 30 -Q 20 -l `"$posPath`" -f `"$refPath`" $bamlist_space | `"$SED`" `"s/chr//`""
    }
    $suffix = " | `"$GAWK`" `"$gawkExpr`" | `"$USORT`" -t `"$tab`" -k1,1n -k2,2n | `"$PILEUPCALLER`" --randomHaploid --sampleNames $array2_string --samplePopName $PopName -f `"$snpPath`" -p `"$outPath`" > `"$statsPath`" 2>&1"

    $batchContent = "@echo off`r`n$prefix$suffix`r`n"
    $tempBatchPath = [System.IO.Path]::GetTempFileName() + ".bat"
    Set-Content -Path $tempBatchPath -Value $batchContent -Encoding ASCII

    try {
        & $tempBatchPath
    }
    finally {
        Remove-Item $tempBatchPath -ErrorAction SilentlyContinue
    }
}

# Function to process files in batches when a single command would exceed cmd's length limit
function Invoke-FilesInBatch {
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

        try {
            Invoke-Pipeline -FileList $batch -ReferenceType $ReferenceType -PopName $PopName -OutputName $OutputName -BatchSuffix "_batch$currentBatch"
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
    $response = Read-Host "Do you want to download it? (Y/N)"
    if ($response -eq 'Y') {
        Write-Host "Downloading"
        Push-Location reference
        & $CURL --insecure -O "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
        Write-Host "Decompressing"
        & $BGZIP -d -@2 hs37d5.fa.gz
        Write-Host "Indexing hs37d5.fa"
        & $SAMTOOLS faidx hs37d5.fa
        Pop-Location
    }
    else {
        exit 0
    }
    Write-Host ""
}

# Check for hg19.fa
if (-not (Test-Path "reference\hg19.fa")) {
    Write-Host "`nYou need to have hg19.fa in the reference directory"
    $response = Read-Host "Do you want to download it? (Y/N)"
    if ($response -eq 'Y') {
        Write-Host "Downloading"
        Push-Location reference
        & $CURL --insecure -O "https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz"
        Write-Host "Decompressing"
        & $BGZIP -d -@2 hg19.fa.gz
        Write-Host "Indexing hg19.fa"
        & $SAMTOOLS faidx hg19.fa
        Pop-Location
    }
    else {
        exit 0
    }
    Write-Host ""
}

# Prompt for Population and Output names
$POPNAME    = Read-Host "Enter Population name"
$OUTPUTNAME = Read-Host "Enter Output name"

# Navigate to source directory and process BAM/CRAM files
Push-Location $Source
$headerCheck = $null
$FLAG        = "NEITHER"
$bamlist1    = @()

Get-ChildItem -Path * -Include "*.bam", "*.cram" | ForEach-Object {
    if ($null -eq $headerCheck) {
        $headerCheck = $_.Name
        $HEADER_TEMP = "header_temp.txt"
        & $SAMTOOLS view -H $headerCheck > $HEADER_TEMP

        if (Select-String -Path $HEADER_TEMP -Pattern (Get-Content (Join-Path $BASEDIR "patternhs37d5.txt"))) {
            $FLAG = "HS37D5"
        }
        elseif (Select-String -Path $HEADER_TEMP -Pattern (Get-Content (Join-Path $BASEDIR "patternhg19.txt"))) {
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

Write-Host "Attempting to process all files at once..."

try {
    Invoke-Pipeline -FileList $bamlist1 -ReferenceType $FLAG -PopName $POPNAME -OutputName $OUTPUTNAME
}
catch {
    Write-Error "Command failed. Error details: $_"
    Write-Host "Attempting batch processing instead..."
    $batchSuccess = Invoke-FilesInBatch -FileList $bamlist1 -ReferenceType $FLAG -PopName $POPNAME -OutputName $OUTPUTNAME
    if (-not $batchSuccess) {
        Write-Error "Both single command and batch processing failed"
        Pop-Location
        exit 1
    }
}

Pop-Location
Write-Host "Done!"
Read-Host "Press Enter to continue"
