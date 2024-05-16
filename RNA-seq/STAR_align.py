# !/usr/bin/env python
#_*_coding:utf-8 _*_
import multiprocessing.pool
import pathlib
import os
import shutil
import argparse
import sys
import subprocess
import glob
import pysam
import multiprocessing
from functools import partial

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

# Module:: Generate bash command
def create_bash(script_name, command):
    with open(script_name, 'w') as file:
        file.write("#!/bin/bash\n")
        file.write(command)
    file.close()
    
# Module:: Execute bash command
def execute_command(script_name):
    status, output = subprocess.getstatusoutput("bash {}".format(script_name))
    if status == 0:
        print("{} ({})".format(output, script_name))
        return output.split()[-1]
    else:
        print("Error when executing the following command: {} {} {}".format(status, output, script_name))


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
          
    create_bash("{}/tmp/fastp_trim_{}.sh".format(output_dir,sample_ID),command)
    return execute_command("{}/tmp/fastp_trim_{}.sh".format(output_dir,sample_ID))

# Module:: STAR alignment
def STAR_align(input_dir, sample_ID, output_dir, gtf, 
               star_path, star_index, paired, samtype):
    
    if paired:
        command = "{} --runMode alignReads ".format(star_path) 
        command += "--runThreadN 24 --genomeDir {} ".format(star_index)
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
        command += "--runThreadN 24 --genomeDir {} ".format(star_index)
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
        
    create_bash("{}/tmp/STAR_align_{}.sh".format(output_dir,sample_ID),command)
    return execute_command("{}/tmp/STAR_align_{}.sh".format(output_dir,sample_ID))

# Module:: Sort and index BAM files
def sort_index_bam(output_dir, sample_ID, SAMtype):
    
    if SAMtype == 'BAM Unsorted':
        sam_dir = pathlib.Path(output_dir) / 'bam_unsorted'
        sample_dir = pathlib.Path(sam_dir) / sample_ID
        bam_file = pathlib.Path(sample_dir) / '{}Aligned.out.bam'.format(sample_ID)
        bam_file_sort = bam_file.with_suffix(".sorted.bam")
        bai_file = bam_file_sort.with_suffix(".bai")
        print("Need sorting first! Perform sorting using samtools......\n")
        
        command = "samtools sort -@ 16 -m 2G -o {} {} && ".format(bam_file_sort, bam_file)
        command += "samtools index -b -@ 16 {} {}".format(bam_file_sort, bai_file)
        
    elif SAMtype == 'SAM':
        sam_dir = pathlib.Path(output_dir) / 'sam'
        sample_dir = pathlib.Path(sam_dir) / sample_ID
        bam_file = pathlib.Path(sample_dir) / '{}Aligned.out.sam'.format(sample_ID)
        bam_file_sort = bam_file.with_suffix(".sorted.sam")
        bai_file = bam_file_sort.with_suffix(".bai")
        print("Need sorting first! Perform sorting using samtools......\n")
        
        command = "samtools sort -@ 16 -m 2G -o {} {} && ".format(bam_file_sort, bam_file)
        command += "samtools index -b -@ 16 {} {}".format(bam_file_sort, bai_file)
             
    elif SAMtype == 'BAM SortedByCoordinate':
        sam_dir = pathlib.Path(output_dir) / 'bam_sorted'
        sample_dir = pathlib.Path(sam_dir) / sample_ID
        bam_file = pathlib.Path(sample_dir) / '{}Aligned.sortedByCoord.out.bam'.format(sample_ID)
        bai_file = pathlib.Path(bam_file.with_suffix(".bai"))
        print("BAM has been sorted by STAR. Perform indexing using samtools......\n")
        command = "samtools index -b -@ 16 {} {} ".format(bam_file, bai_file)
        

    create_bash("{}/tmp/sort_index_{}.sh".format(output_dir,sample_ID),command)
    return execute_command("{}/tmp/sort_index_{}.sh".format(output_dir,sample_ID))


# Module:: Process each sample
def process_sample(input_dir, sample_ID, output_dir, gtf, 
                   star_path, star_index, paired, SAMtype):
        
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
    #print(sam_dir)
    sample_dir = pathlib.Path(sam_dir) / sample_ID
    sample_dir.mkdir(parents=True, exist_ok=True)
    #print(sample_dir)
    
    # Check if trimming is performed
    if paired:
        check_trim = glob.glob('{}/{}_1_trimmed.fastq.gz'.format(input_dir,sample_ID))
    else:
        check_trim = glob.glob('{}/{}_trimmed.fastq.gz'.format(input_dir,sample_ID))
    
    # Check if alignment is performed
    check_bam = glob.glob('{}/*.bam'.format(sample_dir))
    '''
    # Check if sort is performed
    if SAMtype == 'BAM Unsorted':
        check_sorted_bam = glob.glob('{}/*.sorted.bam'.format(sample_dir))
    
    elif SAMtype == 'BAM SortedByCoordinate':
        check_sorted_bam = check_bam
    
    elif SAMtype == 'SAM':
        check_sorted_bam = glob.glob('{}/*.sorted.sam'.format(sample_dir))
    '''
    #print((not check_trim), (not check_bam))
    if not check_trim:
        print("Performing fastp trimming and STAR alignment for {}......\n".format(sample_ID))
        fastp_trim(input_dir, paired, sample_ID, output_dir)
        STAR_align(input_dir, sample_ID, output_dir, gtf, star_path, star_index, paired, SAMtype)
        #sort_index_bam(output_dir, sample_ID, SAMtype)
        
    elif not check_bam:
        print("Trimming has been performed. Performing STAR alignment for {}......\n".format(sample_ID)) 
        STAR_align(input_dir, sample_ID, output_dir, gtf, star_path, star_index, paired, SAMtype)
        #sort_index_bam(output_dir, sample_ID, SAMtype)
    '''            
    elif not check_sorted_bam:
        print("Fastp trimming and STAR alignment have been performed for {}.\n".format(sample_ID))
        sort_index_bam(output_dir, sample_ID, SAMtype)
    
    else:    
        print("Fastp trimming, STAR alignment and Sorting have been performed for {}.\n".format(sample_ID))
    '''
    # After alignment, perform sorting and indexing for SAM/BAM files. 
    sort_index_bam(output_dir, sample_ID, SAMtype)



# Main Module:: RNA-seq pipeline
def rna_seq_pipe(input_dir, output_dir, paired, gtf, star_path, star_index, SAMtype, keep):
    # Get iteratable sample names 
    names = get_sample_ID(input_dir, paired) # Example output: names = ['SRR123', 'SRR234']
    
    # Start multi-processing for fastp trimming and STAR alignment
    pool = multiprocessing.Pool(processes=4)
    for name in names:  
        pool.apply_async(process_sample, args=(input_dir, name, output_dir, gtf, 
                                               star_path, star_index, paired, SAMtype))
    print('Multi-tasks for RNA-seq pipeline in progress...\n')
    pool.close()
    pool.join()
    print('Pipeline completed!')

    # Remove tmp folder if keep is False
    work_dir = pathlib.Path(output_dir) / 'tmp'
    
    if not keep:
        shutil.rmtree(work_dir)

# Execute the main function
rna_seq_pipe(args.input, args.output, args.paired, args.gtf, 
             args.star_path, args.star_index, args.samtype, args.keep)

# End of the script