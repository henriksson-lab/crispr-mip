#https://umi-tools.readthedocs.io/en/latest/API.html
#pip install umi_tools
import gzip
import io
import pandas as pd
import numpy as np
import sys

print("subsamp for pcr")

########################## Get tables of expected gRNA sequences
grnalist = pd.read_csv("/home/mahogny/mystore/dataset/crispr_padlock/lib/brunello.csv",sep="\t")
set_grna_brunello = set(grnalist["sgRNA Target Sequence"])  #20bp

grnalist = pd.read_csv("/home/mahogny/mystore/dataset/crispr_padlock/lib/kinome.csv",sep="\t")
set_grna_kinome = set(grnalist["sgRNA Target Sequence"])  #20bp


all_set_grna = dict({
    "Brunello":set_grna_brunello,
    "Brunello-kinome":set_grna_kinome,
    "GFPg1":set_grna_brunello  #hack
})


########################## Read metadata to figure out what files we have
samplemeta = pd.read_csv("/home/mahogny/mystore/dataset/crispr_padlock/samplemeta.csv",sep="\t")
samplemeta = samplemeta[samplemeta["protocol"]=="PCR"]

rootdir = "/home/mahogny/mystore/dataset/crispr_padlock/fastq/"


#########
# Function to extract the sgRNA, which is in front of afterseq
# 20bp in front of "GTTTTAGAGCTAGAAA"
afterseq = "GTTTTAGAGC"
def getgrna(oneread):
    p1 = oneread.split(afterseq)[0]
    return p1[(len(p1)-20):len(p1)]    
#grna = [getgrna(r) for r in allreadsR1]


def readBufferToGrnaOLD(seqR1):
    with io.TextIOWrapper(io.BufferedReader(gzip.open(seqR1))) as buffer:
        lines=buffer.readlines()
    allreads = [getgrna(lines[x+1]) for x in range(0, len(lines), 4)]
    print("Read "+str(len(allreads))+" reads")
    return allreads

def readBufferToGrnaNEW(seqR1):
    with io.TextIOWrapper(io.BufferedReader(gzip.open(seqR1))) as buffer:
        lines=buffer.readlines()
    allreads = [getgrna(lines[x+1]) for x in range(0, len(lines), 4)]
    print("Read "+str(len(allreads))+" reads")
    return allreads


#########
# Read FASTQ files, extract gRNA and UMI from reads.
# Discard reads that do not perfectly match brunello list
def readFastqPCR(seqR1, seqR2, seqR1b, seqR2b, set_grna):


    ### Extract all reads from fastq
    grna = readBufferToGrnaOLD(seqR1)

    ### If resequenced, add additional reads
    if seqR1b!="":
         grna = grna + readBufferToGrnaOLD(seqR1b)  ########## bug. need to return it all with seqR1b instead

    ### Only keep grnas in table
    bad_grna = [x for x in grna if x not in set_grna]
    grna = [x for x in grna if x in set_grna]


    table_reads = pd.DataFrame({"grna": bad_grna})
    table_reads.to_csv("/home/mahogny/mystore/dataset/crispr_padlock/prob.csv", index=False)


    table_reads = pd.DataFrame({"grna": grna})

    print("=== after reading, total grnas # "+str(table_reads.shape[0]))

    return table_reads



subsample_size_list = [0,100,250, 500,750, 1000,2000, 3000, 4000, 5000,7000, 10000, 15000, 20000, 25000, 50000 ] + [round(x*1000000) for x in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 15]]

i = int(sys.argv[1])

## lib PCR is 4..7  and not getting it all??

#################################### Process all input files: PCR
#for i in range(samplemeta.shape[0]):
if True:
    cursamp = list(samplemeta["filename"])[i]
    cursamp_re = list(samplemeta["reseq_sample"])[i]
    print(i)
    print(cursamp)
    print(cursamp_re)

    if type(cursamp_re) == str:
        fR1 = rootdir+cursamp_re
        fR2 = rootdir+cursamp_re.replace("_R1","_R2")
    else:
        fR1 = ""
        fR2 = ""

    print("additional reads")
    print(fR1)
    print(fR2)


    ## Read all the data, extract grna,umi  
    table_reads = readFastqPCR(
        rootdir+cursamp,
        rootdir+cursamp.replace("_R1","_R2"),
        fR1,
        fR2,
        all_set_grna[list(samplemeta["library"])[i]]
    )

    
    
    for cursize in subsample_size_list:
        print("===Attempting size: ", str(cursize))
        if cursize < table_reads.shape[0]:
            num_subsamp = 5
            if cursize==0:
                num_subsamp = 1

            for subsamp_rng in range(0,num_subsamp):
                print("do "+str(cursize)+" # "+str(subsamp_rng))
                np.random.seed(subsamp_rng)

                if cursize==0:
                    sub_table_reads = table_reads
                    cursize = "allreads"
                else:
                    sub_table_reads = table_reads.sample(n=cursize, replace=False, random_state=subsamp_rng)

                counts = sub_table_reads['grna'].value_counts()
                d = pd.DataFrame({"grna":counts.index, "cnt":counts})
                d.to_csv("/home/mahogny/mystore/dataset/crispr_padlock/subsamp/"+cursamp+"#"+str(cursize)+"#"+str(subsamp_rng), index=False)

    
