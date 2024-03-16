# aDNA to dataset

Added a reference downloader to both versions.

You need to have hs37d5.fa, hs37d5.fa.fai, hg19.fa and hg19.fa.fai in reference subdirectory.

hs37d5.fa.gz + fai can be downloaded from:

https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/

hg19.fa.gz + fai can be downloaded from:

https://omic.tech/ftp/public/3dsnp/assembly/

Notice! You have to decompress and index with samtools faidx both of the above references!

Simple instructions:

Copy all of the bams you want to include in source directory and run the program. The result will be in target subdirectory.
