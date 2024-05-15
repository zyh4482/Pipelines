# !/usr/bin/env python
#_*_coding:utf-8 _*_
import pathlib
import os
import shutil
import argparse
import sys
import subprocess
import glob
import pysam
#import multiprocessing

# Define the parser
parser = argparse.ArgumentParser(description='Perform STAR alignment for RNA-seq')
parser.add_argument('--input', '-i', action="store",dest="input", help="Enter Input directory", required=True)
parser.add_argument('--out', '-o', action="store",dest="output", help="Enter Output directory", required=True)
parser.add_argument('--pair', '-p', action="store_true",dest="paired", 
                    help="Paired reads or not", required=False) # If this flag is not set, it will be False. Reads will be taken as single-end reads
parser.add_argument('--gtf', '-g', action="store",dest="gtf", help="Enter GTF directory", required=True)
parser.add_argument('--star_path', '-s', action="store",dest="star_path", help="Enter STAR executable directory", required=True)
parser.add_argument('--star_index', '-I', action="store",dest="star_index", help="Enter STAR Index directory", required=True)
parser.add_argument('--SAMtype', '-t', action="store",dest="samtype", 
                    choices = ['BAM Unsorted', 'BAM SortedByCoordinate', 'SAM'],
                    help="Decide which type of aligned file as output", required=True)
parser.add_argument('--keep', '-k', action="store_true",dest="keep", 
                    help="Keep tmp scripts in tmp folder. Default: False", required=False) # If this flag is not set, it will be False. tmp folder will be removed

args = parser.parse_args()

# Automatically print help message if argument is not included
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

#pathlib.Path(args.output+'/tmp').mkdir(parents=True, exist_ok=True)


# Module:: Generate bash command
def create_bash(file_name, command):
    with open(file_name, 'w') as file:
        file.write("#!/bin/bash\n")
        file.write(command)
    file.close()
    
'''
# Module:: Execute bash command
def submit_job(file_name):
    status, output = subprocess.getstatusoutput("bash {}".format(file_name))
    if status == 0:
        print("{} ({})".format(output, file_name))
        return output.split()[-1]
    else:
        print("Error when executing the following command: {} {} {}".format(status, output, file_name))
'''

# Module:: Get sample ID of Fastq files
def get_sample_ID(input_dir, paired):
    if paired:
        names = glob.glob(input_dir + '/*_1.f*')
        sample_ID = [pathlib.Path(name).name.split('_1')[0] for name in names]
        return sample_ID
    else: 
        names = glob.glob(input_dir + '/*.f*')
        sample_ID = [pathlib.Path(name).name.split('.f')[0] for name in names]
        return sample_ID

# Module:: fastp trimming
def fastp_trim(input_dir, paired, sample_ID,output_dir):
    suffixes = ['.fastq.gz', '.fq.gz', '.fq', '.fastq']
    
    if paired: 
        
        for suffix in suffixes:
            file1 = os.path.join(input_dir, "{}_1{}".format(sample_ID, suffix)) 
            # example:file1 = "/home/user/data/sample_1.fastq.gz"
            file2 = os.path.join(input_dir, "{}_2{}".format(sample_ID, suffix))
            # example: file2 = "/home/user/data/sample_2.fastq.gz"
            if os.path.exists(file1) and os.path.exists(file2):
                command = "fastp --thread 16 --qualified_quality_phred 20 -f 15 --cut_right --detect_adapter_for_pe "
                command += "--cut_window_size 4 --cut_mean_quality 20 --unqualified_percent_limit 40 --n_base_limit 5 "
                command += "--compression 4 "
                command += "-i {} -I {} ".format(file1, file2)       
                command += "-o {}/{}_1_trimmed.fastq.gz -O {}/{}_2_trimmed.fastq.gz ".format(input_dir, sample_ID, input_dir, sample_ID)   
                command += "--json {}/{}_fastp.json --html {}/{}_fastp.html\n".format(input_dir, sample_ID, input_dir, sample_ID)

            else:
                #print("Error: {} or {} does not exist".format(file1, file2)) # No need to print this
                continue

    else: 
        
        for suffix in suffixes:
            file = os.path.join(input_dir, "{}{}".format(sample_ID, suffix)) 
            # example:file1 = "/home/user/data/sample.fastq.gz"
            if os.path.exists(file):
                command = "fastp --thread 16 --qualified_quality_phred 20 -f 15 --cut_right --detect_adapter_for_pe "
                command += "--cut_window_size 4 --cut_mean_quality 20 --unqualified_percent_limit 40 --n_base_limit 5 "
                command += "--compression 4 "
                command += "-i {} ".format(file)       
                command += "-o {}/{}_trimmed.fastq.gz ".format(input_dir, sample_ID)   
                command += "--json {}/{}_fastp.json --html {}/{}_fastp.html\n".format(input_dir, sample_ID, input_dir, sample_ID)
                
            else:
                #print("Error: {} does not exist".format(file)) # No need to print this
                continue
          
    return create_bash("{}/tmp/fastp_trim_{}.sh".format(output_dir,sample_ID),command)
    #return submit_job("{}/tmp/fastp_trim_{}.sh".format(output_dir,prefix))

# Module:: STAR alignment
def STAR_align(input_dir, sample_ID, output_dir, gtf, 
               star_path, star_index, paired, samtype):
    
    if paired:
        command = "{} --runMode alignReads ".format(star_path) 
        command += "--runThreadN 48 --genomeDir {} ".format(star_index)
        command += "--readFilesCommand zcat --twopassMode Basic --outSAMtype {} ".format(samtype)
        command += "--sjdbGTFfile {} ".format(gtf)
        command += "--outSAMstrandField intronMotif " # Important XS tag for alternative splicing analysis
        command += "--outSAMattributes All "
        command += "--readFilesIn {}/{}_1_trimmed.fastq.gz {}/{}_2_trimmed.fastq.gz ".format(input_dir, sample_ID, input_dir, sample_ID)
        
        if samtype == 'BAM Unsorted':
            command += "--outFileNamePrefix {}/bam_unsorted/{}/{} ".format(output_dir, sample_ID, sample_ID)
            
        elif samtype == 'BAM SortedByCoordinate':
            command += "--outFileNamePrefix {}/bam_sorted/{}/{} ".format(output_dir, sample_ID, sample_ID)
            
        elif samtype == 'SAM':
            command += "--outFileNamePrefix {}/sam/{}/{} ".format(output_dir, sample_ID, sample_ID)
        
        command += "--outSAMattrRGline ID:{} SM:{} LB:{} PL:ILLUMINA ".format(sample_ID, sample_ID, sample_ID)

    else:
        command = "{} --runMode alignReads ".format(star_path) 
        command += "--runThreadN 48 --genomeDir {} ".format(star_index)
        command += "--readFilesCommand zcat --twopassMode Basic --outSAMtype {} ".format(samtype)
        command += "--sjdbGTFfile {} ".format(gtf)
        command += "--outSAMstrandField intronMotif " # Important XS tag for alternative splicing analysis
        command += "--outSAMattributes All "
        command += "--readFilesIn {}/{}_trimmed.fastq.gz ".format(input_dir, sample_ID)
        
        if samtype == 'BAM Unsorted':
            command += "--outFileNamePrefix {}/bam_unsorted/{}/{} ".format(output_dir, sample_ID, sample_ID)
            
        elif samtype == 'BAM SortedByCoordinate':
            command += "--outFileNamePrefix {}/bam_sorted/{}/{} ".format(output_dir, sample_ID, sample_ID)
            
        elif samtype == 'SAM':
            command += "--outFileNamePrefix {}/sam/{}/{} ".format(output_dir, sample_ID, sample_ID)
        
        command += "--outSAMattrRGline ID:{} SM:{} LB:{} PL:ILLUMINA ".format(sample_ID, sample_ID, sample_ID)

    return create_bash("{}/tmp/STAR_align_{}.sh".format(output_dir,sample_ID),command)
    #return submit_job("{}/tmp/STAR_align_{}.sh".format(output_dir,sample_ID))


# Module:: Main function
def main(input_dir, output_dir, paired, gtf, star_path, star_index, SAMtype, keep):
    
    # Create directory to store verbose tmp scripts
    work_dir = pathlib.Path(output_dir) / 'tmp'
    work_dir.mkdir(parents=True, exist_ok=True)
    
    # Create directory to store SAM/BAM files
    if SAMtype == 'BAM Unsorted':
        sam_dir = pathlib.Path(output_dir) / 'bam_unsorted'
        sam_dir.mkdir(parents=True, exist_ok=True)
    
    elif SAMtype == 'BAM SortedByCoordinate':
        sam_dir = pathlib.Path(output_dir) / 'bam_sorted'
        sam_dir.mkdir(parents=True, exist_ok=True)
    
    elif SAMtype == 'SAM':
        sam_dir = pathlib.Path(output_dir) / 'sam'
        sam_dir.mkdir(parents=True, exist_ok=True)
    
    #star_path = 'STAR' # Manually define the path to STAR # If install with conda, just type-in STAR
    #star_index = '/home/zyh/ref/star_index/GRCh38'
    #gtf = '/home/zyh/ref/gencode.v40.annotation.gtf'
    
    names = get_sample_ID(input_dir, paired) # Example output: names = ['SRR123', 'SRR234']
    
    # For each sample, create shell scripts for fastp trimming and STAR alignment
    # Shell script is stored in the tmp folder under the given output directory
    
    for sample_ID in names:
        
        sample_dir = pathlib.Path(sam_dir) / sample_ID
        sample_dir.mkdir(parents=True, exist_ok=True)
        
        if paired:
            check_trim = glob.glob('{}/{}_1_trimmed.fastq.gz'.format(input_dir,sample_ID))
        else:
            check_trim = glob.glob('{}/{}_trimmed.fastq.gz'.format(input_dir,sample_ID))
            
        check_bam = glob.glob('{}/*.bam'.format(sample_dir))
        
        if not check_trim:
            print("Performing fastp trimming and STAR alignment for {}......\n".format(sample_ID))
            fastp_trim(input_dir, paired, sample_ID, output_dir)
            STAR_align(input_dir, sample_ID, output_dir, gtf, star_path, star_index, paired, SAMtype)
        
        elif not check_bam:
            print("Trimming has been performed.Performing STAR alignment for {}......\n".format(sample_ID)) 
            STAR_align(input_dir, sample_ID, output_dir, gtf, star_path, star_index, paired, SAMtype)
        
        else:
            print("Fastp trimming and STAR alignment have been performed for {}.\n".format(sample_ID))
            continue   

        '''
        if SAMtype == 'BAM Unsorted':
            bam_file = pathlib.Path(sample_dir) / '{}Aligned.out.bam'.format(sample_ID)
            print("Need sorting first! Perform sorting using samtools......\n")
            
            pysam.sort("-@" "16", "-m", "2G","-o", bam_file.replace(".bam", ".sorted.bam"), bam_file)
            
            bam_file_sort = bam_file.replace(".bam", ".sorted.bam")
            pysam.index(bam_file_sort, bam_file_sort.replace(".sorted.bam", ".sorted.bai"),"-b","-@","16")
            
        elif SAMtype == 'SAM':
            bam_file = pathlib.Path(sample_dir) / '{}Aligned.out.sam'.format(sample_ID)
            print("Need sorting first! Perform sorting using samtools......\n")
            
            pysam.sort("-@" "16", "-m", "2G","-o", bam_file.replace(".sam", ".sorted.sam"), bam_file)
            
            bam_file_sort = bam_file.replace(".sam", ".sorted.sam")
            pysam.index(bam_file_sort, bam_file_sort.replace(".sorted.sam", ".sorted.bai"),"-b","-@","16")
            
        elif SAMtype == 'BAM SortedByCoordinate':
            bam_file = pathlib.Path(sample_dir) / '{}Aligned.sortedByCoord.out.bam'.format(sample_ID)
            print("BAM has been sorted by STAR. Perform indexing using samtools......\n")
            pysam.index(bam_file, bam_file.replace(".bam", ".bai"),"-b","-@","16")
        '''

    # Find all fastp_trim shell scripts and execute them
    fastp_trim_scripts = glob.glob('{}/fastp_trim*.sh'.format(work_dir))
    procs = [ subprocess.Popen("bash "+i, shell=True) for i in fastp_trim_scripts]
    for p in procs:
        p.wait()
        
    # Find all STAR_align shell scripts and execute them
    STAR_align_scripts = glob.glob('{}/STAR_align*.sh'.format(work_dir))
    #print([i for i in STAR_align_scripts])
    procs = [ subprocess.Popen("bash "+i, shell=True) for i in STAR_align_scripts]
    for p in procs:
        p.wait()
        
    for sample_ID in names:
        if SAMtype == 'BAM Unsorted':
            bam_file = pathlib.Path(sample_dir) / '{}Aligned.out.bam'.format(sample_ID)
            print("Need sorting first! Perform sorting using samtools......\n")
            
            pysam.sort("-@" "16", "-m", "2G","-o", bam_file.replace(".bam", ".sorted.bam"), bam_file)
            
            bam_file_sort = bam_file.replace(".bam", ".sorted.bam")
            pysam.index(bam_file_sort, bam_file_sort.replace(".sorted.bam", ".sorted.bai"),"-b","-@","16")
            
        elif SAMtype == 'SAM':
            bam_file = pathlib.Path(sample_dir) / '{}Aligned.out.sam'.format(sample_ID)
            print("Need sorting first! Perform sorting using samtools......\n")
            
            pysam.sort("-@" "16", "-m", "2G","-o", bam_file.replace(".sam", ".sorted.sam"), bam_file)
            
            bam_file_sort = bam_file.replace(".sam", ".sorted.sam")
            pysam.index(bam_file_sort, bam_file_sort.replace(".sorted.sam", ".sorted.bai"),"-b","-@","16")
            
        elif SAMtype == 'BAM SortedByCoordinate':
            bam_file = pathlib.Path(sample_dir) / '{}Aligned.sortedByCoord.out.bam'.format(sample_ID)
            print("BAM has been sorted by STAR. Perform indexing using samtools......\n")
            pysam.index(bam_file, bam_file.replace(".bam", ".bai"),"-b","-@","16")    
    
    # Remove tmp folder if keep is False
    if not keep:
        shutil.rmtree(work_dir)


# Execute the main function
main(input_dir=args.input, output_dir=args.output, 
     paired=args.paired, gtf=args.gtf, star_path=args.star_path, 
     star_index=args.star_index, SAMtype=args.samtype, keep=args.keep)
# End of the script