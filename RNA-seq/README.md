# RNA-seq

This pipeline is built to perform RNA-seq processing and prepare necessary data for alternative splicing detection under python. 

## Reads alignment ##
`STAR_align.py` will automatically perform trimming, alignment, sorting and indexing for bulk RNA-seq. 

### Dependencies: 

`fastp 0.20.0`


`STAR 2.7.10a`

### Usage:

Example code: 
```
python /home/usr/script/STAR_align.py -i /home/data/fastq \
-o /home/usr/data/bam --pair --star_path STAR --star_index /home/usr/ref/star_index/GRCh38/gencode40 \
--gtf /home/usr/ref/gencode.v40.annotation.gtf --SAMtype 'BAM SortedByCoordinate'
```

#### Arguments
`STAR_align.py` takes following arguments:

`--input`,`-i`: Enter input directory. This directory is the one store your raw `fastq` or `fastq.gz` files. 

`--out`, `-o`: Enter output directory. 

`--pair`, `-p`: If this flag is added, it will recognize files as single end sequencing data. 

`-gtf`, `-g`: Enter path to GTF. e.g. `/home/usr/ref/gencode40.gtf`

`--star_path`, `-s`: Enter the executable of STAR. For users who installed STAR using conda, please just enter `STAR` for this argument. 

`--star_index`, `-I`: Enter STAR index directory. 

`--SAMtype`, `-t`: Decide which type of SAM/BAM file you wish to create. Options: `SAM`, `BAM Unsorted`, `BAM SortedByCoordinate`. 

`--keep`, `-k`: If this flag is added, keep the verbose scripts stored in `/tmp` folder. 


### Description

This script will automatically perform trimming, alignment, sorting and indexing for bulk RNA-seq data. 

If your input directory is `/home/usr/data/fastq/`: 

Trimmed `*_trimmed.fastq.gz` files will be stored under the same directory of the raw data. 

Based on different flags you give to `--SAMtype`, it will create a folder beside your input directory named either `bam` or `sam`. 

For example, if we use `--SAMtype BAM Unsorted`, you will find your bam files stored under `/home/usr/data/bam/bam_unsorted/`.
Another folder named `/tmp/` will be created to store temporary shell scripts for trimming and alignment. 

If you choose `SAM` or `BAM Unsorted`, this program will automatically perform sorting and indexing. `*.sorted.bam` and `*.bai` will be stored in the same folder. 

#### Special note: 
You can add flag `--keep` to keep shell scripts stored in `/tmp`. By default, these files will be deleted.
If you did not delete the `/tmp` folder that holds shell scripts of all the samples, this program will fail to detect file existence from `fastp` and `STAR`. 
Therefore, if you rerun the script for some reasons, PLEASE remember to delete the `/tmp` in prior. 

