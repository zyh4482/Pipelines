# !/usr/bin/env python
#_*_coding:utf-8 _*_
import pathlib
import pandas as pd
import subprocess
import argparse
import sys
from subprocess import Popen


# Define the parser
parser = argparse.ArgumentParser(description='Input download directory')
parser.add_argument('--out', '-o', action="store",dest="path", help="enter download directory", required=True)

# Automatically print help message if argument is not included
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)
    
# Now, parse the command line arguments and store the 
# values in the `args` variable
args = parser.parse_args()

#print(args.path)

pathlib.Path(args.path+"/data/fastq").mkdir(parents=True, exist_ok=True) 

data = pd.read_csv("{}/metadata_ena.txt".format(args.path),sep='\t')
data =  data.loc[:,'fastq_aspera'].str.split(';',expand=True)
c1 = "era-fasp@" + data[0]
c2 = "era-fasp@" + data[1]
c1.to_csv("{}/ascp_R1.txt".format(args.path),index=False,header=False,quoting=False)
c2.to_csv("{}/ascp_R2.txt".format(args.path),index=False,header=False,quoting=False)

f = open(pathlib.Path(args.path)/"download_R1.sh",'w')
f.write("# !/bin/bash\n")
f.write("cat "+ args.path +'/ascp_R1.txt'+" | while read id;\n")
f.write("do ascp -QT -I 500m -P33001 -k 1 -i /home/zyh/.aspera/connect/etc/asperaweb_id_dsa.openssh -v ${id} "+args.path+"/data/fastq/;\n")
f.write("done")
f.close()

f = open(pathlib.Path(args.path)/"download_R2.sh",'w')
f.write("# !/bin/bash\n")
f.write("cat " + args.path + '/ascp_R2.txt'+" | while read id;\n")
f.write("do ascp -QT -I 500m -P33001 -k 1 -i /home/zyh/.aspera/connect/etc/asperaweb_id_dsa.openssh -v ${id} "+args.path+"/data/fastq/;\n")
f.write("done")
f.close()

command1 = 'bash '+args.path+'/download_R1.sh'
command2 = 'bash '+args.path+'/download_R2.sh'
commands = [command1,command2]

procs = [ Popen(i, shell=True) for i in commands ]

for p in procs:
   p.wait()