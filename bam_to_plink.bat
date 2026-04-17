@echo off
setlocal EnableDelayedExpansion

REM Capture base directory (where the script is run from)
set "BASEDIR=%CD%"

REM Default values
set "SOURCEDIR=source"

REM Parse command line arguments
:parse_args
if "%~1"=="" goto :args_done
if /i "%~1"=="--source" (
    set "SOURCEDIR=%~2"
    shift
    shift
    goto :parse_args
)
echo Unknown option: %~1
echo Usage: %~n0 [--source path\to\folder]
exit /b 1
:args_done

REM Tool paths
set "SAMTOOLS=%BASEDIR%\winbin\samtools.exe"
set "BGZIP=%BASEDIR%\winbin\bgzip.exe"
set "CURL=%BASEDIR%\winbin\curl.exe"
set "SED=%BASEDIR%\winbin\sed.exe"
set "GAWK=%BASEDIR%\winbin\gawk\gawk.exe"
set "USORT=%BASEDIR%\winbin\usort.exe"
set "PILEUPCALLER=%BASEDIR%\winbin\pileupCaller.exe"

REM Check if hs37d5.fa exists
if not exist reference\hs37d5.fa (
    echo.
    echo You need to have hs37d5.fa in the reference directory
    choice /C YN /M "Do you want to download it"
    if errorlevel 2 exit /b 0
    echo Downloading
    pushd reference
    "%CURL%" --insecure -O https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
    echo Decompressing
    "%BGZIP%" -d -@2 hs37d5.fa.gz
    echo Indexing hs37d5.fa
    "%SAMTOOLS%" faidx hs37d5.fa
    popd
    echo.
)

REM Check if hg19.fa exists
if not exist reference\hg19.fa (
    echo.
    echo You need to have hg19.fa in the reference directory
    choice /C YN /M "Do you want to download it"
    if errorlevel 2 exit /b 0
    echo Downloading
    pushd reference
    "%CURL%" --insecure -O https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz
    echo Decompressing
    "%BGZIP%" -d -@2 hg19.fa.gz
    echo Indexing hg19.fa
    "%SAMTOOLS%" faidx hg19.fa
    popd
    echo.
)

REM Prompt user for Population name
set /p POPNAME=Enter Population name:

REM Prompt user for Output name
set /p OUTPUTNAME=Enter Output name:

REM Navigate to the source directory
pushd "%SOURCEDIR%"

set "headerCheck="
set "bamlist1="
set "FLAG=NEITHER"

for /f "delims=" %%a in ('dir /b *.bam *.cram 2^>nul') do (
    if "!headerCheck!"=="" (
        set "headerCheck=%%a"
        set "HEADER_TEMP=header_temp.txt"
        "%SAMTOOLS%" view -H "!headerCheck!" > "!HEADER_TEMP!"
        for /f "tokens=*" %%i in ('type "!HEADER_TEMP!" ^| findstr /b /g:"%BASEDIR%\patternhs37d5.txt"') do set "FLAG=HS37D5"
        for /f "tokens=*" %%i in ('type "!HEADER_TEMP!" ^| findstr /b /g:"%BASEDIR%\patternhg19.txt"') do set "FLAG=HG19"
        del "!HEADER_TEMP!"
    )
    set "bamlist1=!bamlist1! %%a"
)

REM Build comma-separated sample names (filename without extension)
set "array2="
for %%i in (!bamlist1!) do (
    set "array2=!array2!,%%~ni"
)
set "array2_string=!array2:~1!"

REM Run appropriate pipeline based on reference genome
if "!FLAG!"=="HS37D5" (
    echo Running command
    "%SAMTOOLS%" mpileup -B -q 30 -Q 20 -l "%BASEDIR%\positions\v42.4.1240K.pos" -f "%BASEDIR%\reference\hs37d5.fa" !bamlist1! | "%GAWK%" "BEGIN{OFS=\"\t\"} {if($1==\"X\")$1=23; else if($1==\"Y\")$1=24; else if($1==\"MT\")$1=90; print}" | "%USORT%" -t "	" -k1,1n -k2,2n | "%PILEUPCALLER%" --randomHaploid --sampleNames !array2_string! --samplePopName !POPNAME! -f "%BASEDIR%\positions\v42.4.1240K.snp" -p "%BASEDIR%\target\!OUTPUTNAME!" > "%BASEDIR%\target\!OUTPUTNAME!.stats.txt" 2>&1
) else if "!FLAG!"=="HG19" (
    echo Running command
    "%SAMTOOLS%" mpileup -B -q 30 -Q 20 -l "%BASEDIR%\positions\v42.4.hg19.pos" -f "%BASEDIR%\reference\hg19.fa" !bamlist1! | "%SED%" "s/chr//" | "%GAWK%" "BEGIN{OFS=\"\t\"} {if($1==\"X\")$1=23; else if($1==\"Y\")$1=24; else if($1==\"MT\")$1=90; print}" | "%USORT%" -t "	" -k1,1n -k2,2n | "%PILEUPCALLER%" --randomHaploid --sampleNames !array2_string! --samplePopName !POPNAME! -f "%BASEDIR%\positions\v42.4.1240K.snp" -p "%BASEDIR%\target\!OUTPUTNAME!" > "%BASEDIR%\target\!OUTPUTNAME!.stats.txt" 2>&1
) else (
    echo Reference is not compatible
    popd
    pause
    exit /b 1
)

popd
echo Done!
pause
endlocal
