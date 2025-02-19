#https://umi-tools.readthedocs.io/en/latest/API.html
#pip install umi_tools
import gzip
import io
import pandas as pd
import numpy as np
import sys


print("subsamp for mip")

########################## Get tables of expected gRNA sequences
grnalist = pd.read_csv("/cfs/klemming/home/j/johenr/mystore/crispr_padlock/lib/brunello.csv",sep="\t")  #### wrong!!!
set_grna_brunello = set(grnalist["sgRNA Target Sequence"])  #20bp

grnalist = pd.read_csv("/cfs/klemming/home/j/johenr/mystore/crispr_padlock/lib/kinome.csv",sep="\t")
set_grna_kinome = set(grnalist["sgRNA Target Sequence"])  #20bp

grnalist = pd.read_csv("/cfs/klemming/home/j/johenr/mystore/crispr_padlock/lib/gfp.csv",sep="\t")
set_grna_gfp = set(grnalist["sgRNA Target Sequence"])  #20bp


all_set_grna = dict({
    "Brunello":set_grna_brunello,
    "Brunello-kinome":set_grna_kinome,
    "GFPg1":set_grna_gfp
})


########################## Read metadata to figure out what files we have
samplemeta = pd.read_csv("/cfs/klemming/home/j/johenr/mystore/crispr_padlock/samplemeta.csv",sep="\t")
samplemeta = samplemeta[samplemeta["protocol"]=="padlock"]

rootdir = "/cfs/klemming/home/j/johenr/mystore/crispr_padlock/fastq/"




#########
# Read FASTQ files, extract gRNA and UMI from reads.
# Discard reads that do not perfectly match brunello list
def readFastqMIP(seqR1, seqR2, seqR1b, seqR2b,  set_grna):

    ### Extract all reads from fastq
    with io.TextIOWrapper(io.BufferedReader(gzip.open(seqR1))) as buffer:
        lines=buffer.readlines()
    allreadsR1 = [lines[x+1] for x in range(0, len(lines), 4)]

    with io.TextIOWrapper(io.BufferedReader(gzip.open(seqR2))) as buffer:
        lines=buffer.readlines()
    allreadsR2 = [lines[x+1] for x in range(0, len(lines), 4)]

    ### If resequenced, add additional reads
    if seqR1b!="":
        with io.TextIOWrapper(io.BufferedReader(gzip.open(seqR1))) as buffer:
            lines=buffer.readlines()
        allreadsR1b = [lines[x+1] for x in range(0, len(lines), 4)]

        with io.TextIOWrapper(io.BufferedReader(gzip.open(seqR2))) as buffer:
            lines=buffer.readlines()
        allreadsR2b = [lines[x+1] for x in range(0, len(lines), 4)]

        allreadsR1 = allreadsR1 + allreadsR1b
        allreadsR2 = allreadsR2 + allreadsR2b



    ### Extract the sgRNA, which is the first n letters
    grna = [r[20:(20+20)] for r in allreadsR1]

    ### Extract the UMI, which is the first n letters
    umis = [r[0:13] for r in allreadsR2]

    ## Put read & UMI together. Only keep sgRNAs in library (others missequenced etc).
    ## BWA or other aligner would do this better
    correct_readpairs = [x for x in zip(grna, umis) if x[0] in set_grna]

    table_reads = pd.DataFrame(correct_readpairs,columns=['grna','umi'])

    return table_reads




#########
# For a group of UMIs, look up how many counts there are for each, and return a list of counts.
# We may later wish to omit 1-count UMIs (failed to group with any)
def sum_umicounts(umicount_bytes, onegroup):
    return sum([umicount_bytes[x] for x in onegroup])


#########
# This deduplicates UMIs for a single gRNA; meant to be applied to a data frame with "umi" column
def umicount_per_grna(tab):
    #print(tab)

    ### Count occurence of each UMI according to raw sequencing data
    from collections import Counter
    umicount = dict(Counter(tab["umi"]))
    
    ### UMI-tools expects umis to be byte strings, not regular strings
    umicount_bytes = dict(zip(
        [bytes(k, 'ascii') for k in umicount.keys()],
        umicount.values()))

    ### Gets groups of umis belonging together, but not how many of each umi
    from umi_tools import UMIClusterer
    clusterer = UMIClusterer(cluster_method="directional")
    clustered_umis = clusterer(umicount_bytes, threshold=2) ########### used in the paper. seems safer
    #clustered_umis = clusterer(umicount_bytes, threshold=1) #testing if better for large counts; not good if poor saturation of library
    
    ### Gather counts per umi group
    return [sum_umicounts(umicount_bytes, onegroup) for onegroup in clustered_umis]


#########
# This deduplicates UMIs for all gRNAs in a data frame of (grna,umi). Need not be sorted
def dedup_full_table(table_reads):

    ###### Get a data frame with grna as rownames, and list of counts per umi as one column
    todedup = table_reads.groupby(['grna'])
    deduped = todedup.apply(umicount_per_grna, include_groups=False)    #updated, should be ok
    #DeprecationWarning: DataFrameGroupBy.apply operated on the grouping columns. This behavior is deprecated, and in a future version of pandas the grouping columns will be excluded from the operation. Either pass `include_groups=False` to exclude the groupings or explicitly select the grouping columns after groupby to silence this warning.

    ###### Filter out single-count UMIs, and give a summary count of deduped gRNAs
    def countumi_and_filter(onegroup):
        return len([x for x in onegroup if x>0])

    counttable_deduped = pd.DataFrame({
        "grna":deduped.index, 
        "cnt":[countumi_and_filter(x) for x in deduped]}
    )
    return counttable_deduped




subsample_size_list = [0,100,250, 500,750, 1000,2000, 3000, 4000, 5000,7000, 10000, 15000, 20000, 25000, 50000 ] + [round(x*1000000) for x in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 15]]
#subsample_size_list = [0]

i = int(sys.argv[1])

#################################### Process all input files: MIP
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
    print("Reading all the data")
    table_reads = readFastqMIP(
        rootdir+cursamp,
        rootdir+cursamp.replace("_R1","_R2"),
        fR1,
        fR2,
        all_set_grna[list(samplemeta["library"])[i]]
    )
    print("done reading. num reads in total: "+str(table_reads.shape[0]))
    
    for cursize in subsample_size_list:
        print(cursize)

        if cursize < table_reads.shape[0]:
            print("doing one size")

            num_subsamp = 5
            if cursize==0:
                num_subsamp = 1

            for subsamp_rng in range(0,num_subsamp):
                print("doing one sampling")
                print(str(cursize)+" # "+str(subsamp_rng))
                np.random.seed(subsamp_rng)

                if cursize==0:
                    sub_table_reads = table_reads
                    cursize = "allreads"
                else:
                    sub_table_reads = table_reads.sample(n=cursize, replace=False, random_state=subsamp_rng)

                print("storing non-dedup")
                counts = sub_table_reads[["grna"]].value_counts()
                d = pd.DataFrame({"grna":[x[0] for x in list(counts.index)], "cnt":list(counts)})
                d.to_csv("/cfs/klemming/home/j/johenr/mystore/crispr_padlock/subsamp/"+cursamp+
                         "_nodedup#"+str(cursize)+"#"+str(subsamp_rng), index=False)

                
                if sub_table_reads.shape[0]>0:
                    print("storing dedup")
                    d = dedup_full_table(sub_table_reads) 
                    d.to_csv("/cfs/klemming/home/j/johenr/mystore/crispr_padlock/subsamp/"+cursamp+
                             "_dedup#"+str(cursize)+"#"+str(subsamp_rng), index=False)
                else:
                    print("nothing to dedup "+str(cursize)+"#"+str(subsamp_rng))
        else:
            print("too small")
    
    


print("done")
