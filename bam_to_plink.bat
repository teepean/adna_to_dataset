@echo off
SETLOCAL EnableDelayedExpansion

if not exist reference\hs37d5.fa (
echo .
echo You need to have hs37d5.fa in the reference directory
echo Do you want to download it?
CHOICE /C YN /M "Y/N"
IF ERRORLEVEL == 2  GOTO END
IF ERRORLEVEL == 1 (
echo Downloading
cd reference
..\winbin\curl.exe --insecure -O https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
echo Decompressing
..\winbin\bgzip -d -@2 hs37d5.fa.gz
echo Indexing hs37d5.fa
..\winbin\samtools faidx hs37d5.fa
cd ..
)
echo .
)

if not exist reference\hg19.fa (
echo .
echo You need to have hg19.fa in the reference directory
echo Do you want to download it?
CHOICE /C YN /M "Y/N"
IF ERRORLEVEL == 2  GOTO END
IF ERRORLEVEL == 1 (
echo Downloading
cd reference
..\winbin\curl.exe --insecure -O https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz
echo Decompressing
..\winbin\bgzip -d -@2 hg19.fa.gz
echo Indexing hg19.fa
..\winbin\samtools faidx hg19.fa
cd ..
)
echo .
)
:: Prompt user for Population name
set /p POPNAME=Enter Population name: 

:: Prompt user for Output name
set /p OUTPUTNAME=Enter Output name: 

:: Navigate to the source directory and list all .bam files, storing them in bamlist1
cd source

set "headerCheck="

for /f "delims=" %%a in ('dir /b *.bam *.cram') do (
    if "!headerCheck!"=="" (
        set "headerCheck=%%a"
        set HEADER_TEMP=header_temp.txt
        :: Extract header to temp file
        ..\winbin\samtools view -H !headerCheck! > !HEADER_TEMP!
        :: Initialize flag
        set FLAG=NEITHER

        :: Check patterns
        for /f "tokens=*" %%i in ('type "!HEADER_TEMP!" ^| findstr /b /g:..\patternhs37d5.txt') do set FLAG=HS37D5
        for /f "tokens=*" %%i in ('type "!HEADER_TEMP!" ^| findstr /b /g:..\patternhg19.txt') do set FLAG=HG19
        del "!HEADER_TEMP!"
    )
    set "bamlist1=!bamlist1! %%a"
)

:: Initialize array2 as an empty string
set array2=

:: Loop through bamlist1, remove the .bam extension from each element, and add it to array2
for %%i in (%bamlist1%) do (
    set name=%%~ni
    set array2=!array2!,!name!
)

:: Remove the first comma from array2_string
set array2_string=%array2:~1%

if "!FLAG!"=="HS37D5" (
    echo Running command
    ..\winbin\samtools mpileup -B -q 30 -Q 30 -l ../positions/v42.4.1240K.pos -f ../reference/hs37d5.fa %bamlist1% | ..\winbin\pileupCaller --randomHaploid --sampleNames %array2_string% --samplePopName %POPNAME% -f ../positions/v42.4.1240K.snp -p ../target/%OUTPUTNAME%   > ../target/%OUTPUTNAME%.stats.txt 2>&1
) else if "!FLAG!"=="HG19" (
    echo Running command
    ..\winbin\samtools mpileup -B -q 30 -Q 30 -l ../positions/v42.4.hg19.pos -f ../reference/hg19.fa %bamlist1% | ..\winbin\sed "s/chr//" | ..\winbin\pileupCaller --randomHaploid --sampleNames %array2_string% --samplePopName %POPNAME% -f ../positions/v42.4.1240K.snp -p ../target/%OUTPUTNAME%   > ../target/%OUTPUTNAME%.stats.txt 2>&1
) else (
    echo Reference is not compatible 
    pause
    goto :eof
)

cd ..\target
..\winbin\logdumper.exe
ren console_buffer.txt %OUTPUTNAME%.logfile.txt
cd..

endlocal

pause
