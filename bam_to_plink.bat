@echo off
setlocal enabledelayedexpansion
set PATH=%PATH%;.\winbin

:: Prompt user for Population name
set /p POPNAME=Enter Population name: 

:: Prompt user for Output name
set /p OUTPUTNAME=Enter Output name: 

:: Navigate to the source directory and list all .bam files, storing them in bamlist_space
cd source
set bamlist_space=
for /r %%i in (*.bam) do set bamlist_space=!bamlist_space! "%%i"

:: Initialize array2_string as an empty string
set array2_string=

:: Loop through bamlist_space, remove the .bam extension from each element, and add it to array2_string
for %%i in (!bamlist_space!) do (
    set filename=%%~ni
    set array2_string=!array2_string!,!filename!
)

:: Remove the first comma from array2_string
set array2_string=!array2_string:~1!

set /p choice="Is the reference (hs37d5 or hg19): "
if "%choice%"=="hs37d5" (
    echo Running command
    samtools mpileup -B -q 30 -Q 30 -l ..\positions\v42.4.1240K.pos -f ..\reference\hs37d5.fa !bamlist_space! | pileupCaller --randomHaploid --sampleNames !array2_string! --samplePopName !POPNAME! -f ..\positions\v42.4.1240K.snp -p ..\target\!OUTPUTNAME!
) else if "%choice%"=="hg19" (
    echo Running command
    samtools mpileup -B -q 30 -Q 30 -l ..\positions\v42.4.hg19.pos -f ..\reference\hg19.fa !bamlist_space! | sed "s/chr//" | pileupCaller --randomHaploid --sampleNames !array2_string! --samplePopName !POPNAME! -f ..\positions\v42.4.1240K.snp -p ..\target\!OUTPUTNAME!
) else (
    echo Invalid choice
)
cd ..
