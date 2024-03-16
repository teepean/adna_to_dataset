@echo off
SETLOCAL EnableDelayedExpansion

:: Prompt user for Population name
set /p POPNAME=Enter Population name: 

:: Prompt user for Output name
set /p OUTPUTNAME=Enter Output name: 

:: Navigate to the source directory and list all .bam files, storing them in bamlist1
cd source
for /f "delims=" %%a in ('dir /b *.bam') do (
    set bamlist1=!bamlist1! %%a
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

:: Prompt user for reference choice
set /p choice=Enter your choice (hs37d5 or hg19): 

if "%choice%"=="hs37d5" (
    echo Running command
    ..\winbin\samtools mpileup -B -q 30 -Q 30 -l ../positions/v42.4.1240K.pos -f ../reference/hs37d5.fa.gz %bamlist1% | ..\winbin\pileupCaller --randomHaploid --sampleNames %array2_string% --samplePopName %POPNAME% -f ../positions/v42.4.1240K.snp -p ../target/%OUTPUTNAME%
) else if "%choice%"=="hg19" (
    echo Running command
    ..\winbin\samtools mpileup -B -q 30 -Q 30 -l ../positions/v42.4.hg19.pos -f ../reference/hg19.fa.gz %bamlist1% | ..\winbin\sed "s/chr//" | ..\winbin\pileupCaller --randomHaploid --sampleNames %array2_string% --samplePopName %POPNAME% -f ../positions/v42.4.1240K.snp -p ../target/%OUTPUTNAME%
) else (
    echo Invalid choice
)

cd ..
