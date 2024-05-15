# RNA-seq

This pipeline is built to perform RNA-seq processing and prepare necessary data for alternative splicing detection under python. 

## STAR alignment ##
`STAR_align.py` will automatically perform trimming and alignment for bulk RNA-seq. 

### Dependencies: 

`fastp`

`STAR`

### Usage:

`STAR_align.py` has the following arguments:

`--input`,`-i`: Enter input directory. This directory is the one store your raw `fastq` or `fastq.gz` files. 

`--out`, `-o`: Enter output directory. 

`--pair`, `-p`: If this flag is added, it will recognize files as single end sequencing data. 

`-gtf`, `-g`: Enter path to GTF. e.g. `/home/usr/ref/gencode40.gtf`

`--star_path`, `-s`: Enter the executable of STAR. For users who installed STAR using conda, please just enter `STAR` for this argument. 

`--star_index`, `-I`: Enter STAR index directory. 

`--SAMtype`, `-t`: Decide which type of SAM/BAM file you wish to create. Options: `SAM`, `BAM Unsorted`, `BAM SortedByCoordinate`. 

`--keep`, `-k`: If this flag is added, keep the verbose scripts stored in `/tmp` folder. 


### Description

This script will generate folders as follow:

If your input directory is `/home/usr/data/fastq/`: 

Trimmed `*_trimmed.fastq.gz` files will be stored under the same directory of the raw data. 

Based on different flags you give to `--SAMtype`, it will create a folder beside your input directory named either `bam` or `sam`. 

For example, if we use `--SAMtype BAM Unsorted`, you will find your bam files stored under `/home/usr/data/bam/bam_unsorted/`.
Another folder named `/tmp/` will be created to store temporary shell scripts for trimming and alignment. You can use `--keep` to keep those files. 
By default, these files will be deleted. 

After alignment, you may need to perform sorting and indexing for SAM/BAM. Run `sort_index.py` as below.

--- 
## Sort and index ##


