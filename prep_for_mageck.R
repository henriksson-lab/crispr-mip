
#
# This script generates input files to run mageck
#

library(stringr)
library(dplyr)
library(ggplot2)

################################################################################
###################### Prepare sample meta #####################################
################################################################################


#Input dir
countdir <- "/corgi/otherdataset/crispr_padlock/subsamp/"
#Output dir
dir_mageck <- "/corgi/otherdataset/crispr_padlock/for_mageck"

#What samples we got
samplemeta_all <- read.csv("/husky/martin/CRISPR_padlock_all_fastq/samplemeta.csv",sep="\t")

## Remove samples merged into main sample; not needed?
nrow(samplemeta_all)
reseq_samples <- setdiff(unique(samplemeta_all$reseq_sample),"")
samplemeta_all <- samplemeta_all[!(samplemeta_all$filename %in% reseq_samples),]
nrow(samplemeta_all)

## Split by protocol
samplemeta_pcr   <- samplemeta_all[samplemeta_all$protocol=="PCR",]
samplemeta_padlock_dedup <- samplemeta_all[samplemeta_all$protocol=="padlock",]
samplemeta_padlock_nodedup <- samplemeta_all[samplemeta_all$protocol=="padlock",]


samplemeta_pcr$filename <- paste0(samplemeta_pcr$filename, ".csv")
samplemeta_padlock_dedup$filename <- paste0(samplemeta_padlock_dedup$filename, "_dedup.csv")
samplemeta_padlock_nodedup$filename <- paste0(samplemeta_padlock_nodedup$filename, "_nodedup.csv")

samplemeta_pcr$dedup <- FALSE
samplemeta_padlock_dedup$dedup <- TRUE
samplemeta_padlock_nodedup$dedup <- FALSE

samplemeta <- rbind(
  samplemeta_pcr,
  samplemeta_padlock_dedup,
  samplemeta_padlock_nodedup
)

## Remove bad sample
samplemeta <- samplemeta[!str_detect(samplemeta$filename,"P29968_1005_S5_L001_R1_001"),]




##### Gather sample names for each condition
list_3d_pcr_normoxia <- samplemeta[samplemeta$treatment=="normoxia" & samplemeta$time=="3d" & samplemeta$protocol=="PCR",]$filename
list_3d_mip_dedup_normoxia <- samplemeta[samplemeta$treatment=="normoxia" & samplemeta$time=="3d" & samplemeta$protocol=="padlock" & samplemeta$dedup,]$filename
list_3d_mip_nodedup_normoxia <- samplemeta[samplemeta$treatment=="normoxia" & samplemeta$time=="3d" & samplemeta$protocol=="padlock" & !samplemeta$dedup,]$filename

list_14d_pcr_normoxia <- samplemeta[samplemeta$treatment=="normoxia" & samplemeta$time=="14d" & samplemeta$protocol=="PCR",]$filename
list_14d_mip_dedup_normoxia <- samplemeta[samplemeta$treatment=="normoxia" & samplemeta$time=="14d" & samplemeta$protocol=="padlock" & samplemeta$dedup,]$filename
list_14d_mip_nodedup_normoxia <- samplemeta[samplemeta$treatment=="normoxia" & samplemeta$time=="14d" & samplemeta$protocol=="padlock" & !samplemeta$dedup,]$filename

list_lib_pcr <- samplemeta[samplemeta$treatment=="plasmid_lib" & samplemeta$library=="Brunello-kinome" & samplemeta$protocol=="PCR",]$filename
list_lib_mip_dedup <- samplemeta[samplemeta$treatment=="plasmid_lib" & samplemeta$library=="Brunello-kinome" & samplemeta$protocol=="padlock" & samplemeta$dedup,]$filename
list_lib_mip_nodedup <- samplemeta[samplemeta$treatment=="plasmid_lib" & samplemeta$library=="Brunello-kinome" & samplemeta$protocol=="padlock" & !samplemeta$dedup,]$filename


#PCR library: #note, 4 of them. no reseq
samplemeta[samplemeta$filename %in% list_lib_pcr,] 

#MIP library: #note, 4 of them. no reseq
samplemeta[samplemeta$filename %in% list_lib_mip_dedup,] 

#14d: there is reseq for PCR! but not mip. 
samplemeta[samplemeta$filename %in% list_14d_pcr_normoxia,] 
samplemeta[samplemeta$filename %in% list_14d_mip_nodedup_normoxia,] 


################################################################################
######### Quality control - check sample similarity using UMAP #################
################################################################################

if(FALSE){
  all_cnt <- NULL
  for(fname in list.files(countdir,pattern = "*allreads*")){
    thefile <- file.path(countdir, fname)
    onecnt <- read.csv(thefile)
    onecnt$filename <- fname
    all_cnt[[fname]] <- onecnt
  }
  all_cnt <- do.call(rbind, all_cnt)
  all_cnt$filename <- paste0(str_split_i(all_cnt$f,fixed("#"),1),".csv")
  rownames(all_cnt) <- NULL
  
  #Make a count matrix, normalize to 1
  all_cnt_mat <- acast(all_cnt, filename ~ grna, value.var = "cnt", fill = 0)
  for(i in 1:nrow(all_cnt_mat)){
    all_cnt_mat[i,] <- all_cnt_mat[i,]/sum(all_cnt_mat[i,])
  }

  #Compare samples
  the_umap <- umap(all_cnt_mat)
  
  samplemeta_umap <- merge(samplemeta, 
                           data.frame(
                             umap_x = the_umap$layout[,1],
                             umap_y = the_umap$layout[,2],
                             filename = rownames(the_umap$layout)
                           ))
  
  # list_lib_pcr
  #TODO why no PCR??
  
  ggplot(samplemeta_umap[samplemeta_umap$protocol=="PCR",], aes(umap_x, umap_y, label=samplename, color=paste(protocol,dedup,library))) + geom_point() + geom_text()
  
  ggplot(samplemeta_umap, aes(umap_x, umap_y, label=samplename, color=paste(protocol,dedup,library))) + geom_point() + geom_text()
  ggplot(samplemeta_umap, aes(umap_x, umap_y, label=samplename, color=paste(treatment,time,library))) + geom_point() + geom_text()  
}


if(FALSE){
  ### How similar are two random samples?
  f1 <- read.csv(file.path(countdir, "P29968_1013_S13_L001_R1_001.fastq.gz#3000000#0"))
  colnames(f1) <- c("grna","f1")
  f2 <- read.csv(file.path(countdir, "P29968_1013_S13_L001_R1_001.fastq.gz#3000000#1"))
  colnames(f2) <- c("grna","f2")
  f12 <- merge(f1,f2, all=TRUE)
  ggplot(f12,aes(f1,f2)) + geom_point()
}


################################################################################
### Supplemental figure: comparison of mageck count with our counting method ###
################################################################################

cnt_py <- read.csv(file.path(countdir, "Brunello_kinome_1_S5_merged_R1_001.fastq.gz#allreads#0"))
colnames(cnt_py) <- c("grna","cnt_py")
cnt_mageck <- read.csv("/corgi/otherdataset/crispr_padlock/compare_mageck.out/sample1.count.txt",sep="\t")#[,c(1,3)]
colnames(cnt_mageck) <- c("grna","gene","cnt_mageck")

compcnt <- merge(cnt_py, cnt_mageck)
compcnt$diff <- abs(log10(compcnt$cnt_py) - log10(compcnt$cnt_mageck))

ggplot(compcnt, aes(cnt_py,cnt_mageck, label=gene)) + 
  geom_point(color="gray") + scale_x_log10() + scale_y_log10() +
  geom_text(data=compcnt[compcnt$diff>0.3,]) + 
  theme_bw() +
  xlab("Our counting method")+
  ylab("MAGeCK count")
ggsave("newout/mageck_count.svg", width = 5, height = 5)

compcnt[compcnt$diff>0.3,]


# grna cnt_py   gene cnt_mageck      diff
# 592  AGTCCACCGTCGTGTTCACG     11 FGFRL1         63 0.7579479
# 1693 GATAATATGTTTCATGGACT  20653  MORN2       8740 0.3734717  *********** much higher count in python !
# 2019 GGAGCGTGGAGTCGTCACTT     69   CIB1        205 0.4729048
# 2032 GGAGTGTGACTACCGTCGTG    422  ADCK4       3260 0.8879051
# 2035 GGATAATATGTTTCATGGAC     50  MORN2       7406 2.1706137
# 2068 GGCACTTAGTAAAGTCAGTG     45   TPK1         19 0.3744589
# 2634 TCTCATAAAGCCATTCAATG      3 PAPSS1          1 0.4771213

#Brunello_kinome_1_S5_merged_R1_001.fastq.gz
#mageck count  --list-seq ../mageck_list_grna.txt    --fastq ../fastq/Brunello_kinome_1_S5_merged_R1_001.fastq.gz  --norm-method none
#read.csv("/corgi/otherdataset/crispr_padlock/compare_mageck.out/sample1.count.txt",sep="\t")




################################################################################
######### Prepare grna info in mageck format - functions #######################
################################################################################

grnainfo <- read.csv("/husky/martin/CRISPR_padlock_all_fastq/lib/kinome.csv",sep="\t")

mageck_grna <- unique(data.frame(
  id=grnainfo$sgRNA.Target.Sequence, 
  grna=grnainfo$sgRNA.Target.Sequence,
  gene=grnainfo$Target.Gene.Symbol  
))
mageck_grna$gene[mageck_grna$gene=="Non-Targeting Control"] <- "NTC"
list_control_grna <- data.frame(
  grnaid=mageck_grna$id[mageck_grna$gene=="NTC"]
)

jacks_grna <- data.frame(
  sgRNA=mageck_grna$grna,
  Gene=mageck_grna$gene
)


######## Gather counts
getCounts <- function(mletable, this_samp){
  all_cnt <- NULL
  for(f in mletable$Samples){
    thefile <- file.path(countdir,paste0(f,"#",this_samp))
    if(!file.exists(thefile)){
      print("Missing file")
      print(thefile)
    }
    
    onecnt <- read.csv(thefile)
    onecnt$f <- f
    all_cnt <- rbind(all_cnt, onecnt)
  }
  all_cnt <- reshape2::acast(all_cnt, grna ~ f, value.var = "cnt", fill = 0)
  all_cnt <- as.data.frame(all_cnt)
  all_cnt$id <- rownames(all_cnt)  #using grna as id
  
  all_cnt <- merge(all_cnt,mageck_grna[,c("gene","id")])
  
  all_cnt <- all_cnt[,c("id","gene",mletable$Samples)] #ensure correct column order
  all_cnt
}




######### Set up comparisons, mageck mle
makeMLEtable <- function(listA, listB){
  mletableA <- data.frame(
    Samples=listA,
    baseline=1, 
    eff=1
  )
  mletableB <- data.frame(
    Samples=listB,
    baseline=1, 
    eff=0
  )
  mletable <- rbind(mletableA,mletableB)
  mletable$Samples <- str_remove_all(mletable$Samples,".csv")
  mletable
}




######### Set up comparisons, JACKS
makeJACKStable <- function(listA, listB){
  mletableA <- data.frame(
    Replicate=listA,
    Sample="a"
  )
  mletableB <- data.frame(
    Replicate=listB,
    Sample="b"
  )
  mletable <- rbind(mletableA,mletableB)
  mletable$Replicate <- str_remove_all(mletable$Replicate,".csv")
  mletable
}





######## Figure out what subsamples are available in common
getSampAvail <- function(mletable){
  subsamp_available <- NULL
  subsamp_cnt <- NULL
  for(f in mletable$Samples){
    print(f)
    allf <- list.files(countdir)
    allf <- allf[str_starts(allf,paste0(f,"\\#"))]
    
    these_subsamp <- str_split_fixed(allf,"#",2)[,2]
    
    if(is.null(subsamp_available)){
      subsamp_available <- these_subsamp
      subsamp_cnt <- these_subsamp
    } else {
      subsamp_available <- intersect(subsamp_available, these_subsamp)
      subsamp_cnt <- c(subsamp_cnt, these_subsamp)
    }
  }
  subsamp_available
}

########## Write count files for mageck
writeMageckOut <- function(comparison_name, listA, listB){
  
  #comparison_name <- "14d_pcr"
  #listA <- list_14d_pcr_normoxia
  #listB <- list_lib_pcr
  
  
  ####### Can run comparison if listA and listB has at least one sample!
  
  
  mletable <- makeMLEtable(listA, listB)
  print(mletable)

  jackstable <- makeJACKStable(listA, listB)
  print(jackstable)

  all_cmd_mageck_mle <- c()
  all_cmd_mageck <- c()
  all_cmd_jacks <- c()
  subsamp_available <- getSampAvail(mletable)
  print("samples available")
  print(subsamp_available)
  for(this_samp in subsamp_available){
    #  this_samp <- "100000#0"
    print(this_samp)
    #print(colSums(cnt[,-(1:2)]))
    
    ######## Write out files
    
    cnt <-getCounts(mletable, this_samp)
    write.table(
      cnt, 
      file.path(dir_mageck, paste0(comparison_name, "#",this_samp,".count")),
      sep="\t",
      row.names = FALSE,
      quote = FALSE
    )
    
    write.table(
      mletable,
      file.path(dir_mageck, paste0(comparison_name, "#",this_samp,".mle")),
      sep="\t",
      row.names = FALSE,
      quote = FALSE
    )
    
    write.table(
      jackstable,
      file.path(dir_mageck, paste0(comparison_name, "#",this_samp,"_repmap.tab")),
      sep="\t",
      row.names = FALSE,
      quote = FALSE
    )
    
    
    write.table(
      list_control_grna,
      file.path(dir_mageck, paste0(comparison_name, "#",this_samp,".control")),
      sep="\t",
      col.names = FALSE,
      row.names = FALSE,
      quote = FALSE
    )
    
    write.table(
      mageck_grna,
      file.path(dir_mageck, paste0(comparison_name, "#",this_samp,".grna")),
      sep="\t",
      col.names = FALSE,
      row.names = FALSE,
      quote = FALSE
    )
    
    treatline <- paste(which(mletable$eff==1)-1, collapse = ',')  #updated, rerun all mageck test
    controlline <- paste(which(mletable$eff==0)-1, collapse = ',')  #updated, rerun all mageck test
    
    writeLines(
      treatline,
      file.path(dir_mageck, paste0(comparison_name, "#",this_samp,".treated"))
    )
    
    one_cmd <- paste(
      "mageck test -k ",
      file.path(dir_mageck, paste0(comparison_name, "#",this_samp,".count")),
      "-t", treatline,
      "-c", controlline,
      "--control-sgrna",
      file.path(dir_mageck, paste0(comparison_name, "#",this_samp,".control")),
      "-n",
      file.path(dir_mageck, paste0(comparison_name, "#",this_samp))
    )  
    all_cmd_mageck <- c(all_cmd_mageck, one_cmd)

    ## mageck mle -k ${base_dir}/3_padlock_UMI_${size}.count_normalized.txt 
    #  -d ../design_matrix_padlock_ess.txt -n mageck_output_files_run_3/3_padlock_${size}_mle_UMI 
    #  --control-sgrna ../non_targeting_controls.txt --norm-method median
    

    one_cmd <- paste(
      "mageck mle -k ",
      file.path(dir_mageck, paste0(comparison_name, "#",this_samp,".count")),
      "-d",
      file.path(dir_mageck, paste0(comparison_name, "#",this_samp,".mle")),
      "--control-sgrna",
      file.path(dir_mageck, paste0(comparison_name, "#",this_samp,".control")),
      "-n",
      file.path(dir_mageck, paste0(comparison_name, "#",this_samp,".mageckmle"))
    )  
    all_cmd_mageck_mle <- c(all_cmd_mageck_mle, one_cmd)
    
    
    one_cmd <- paste(
      "python ~/github/JACKS/jacks/run_JACKS.py",
      file.path(dir_mageck, paste0(comparison_name, "#",this_samp,".count")),
      file.path(dir_mageck, paste0(comparison_name, "#",this_samp,"_repmap.tab")),
      "jacks_grnamapping",
      "--common_ctrl_sample=a",
      paste0("--outprefix=",comparison_name, "#",this_samp,"_jacks"),
      "--ctrl_genes=NTC"
    )  
    all_cmd_jacks <- c(all_cmd_jacks, one_cmd)
  }

  replacenullempty <- function(x) if(is.null(x)) "" else x
  
  all_cmd_mageck_mle <- replacenullempty(all_cmd_mageck_mle)
  all_cmd_mageck <- replacenullempty(all_cmd_mageck)
  all_cmd_jacks <- replacenullempty(all_cmd_jacks)
  
  writeLines(all_cmd_mageck_mle, file.path(dir_mageck, paste0(comparison_name, ".mle.sh")))
  writeLines(all_cmd_mageck, file.path(dir_mageck, paste0(comparison_name, ".mageck.sh"))) #mageck test
  writeLines(all_cmd_jacks, file.path(dir_mageck, paste0(comparison_name, ".jacks.sh")))
}



################################################################################
######### Prepare grna info in mageck format - call all function ###############
################################################################################


write.table(
  jacks_grna,
  file.path(dir_mageck, "jacks_grnamapping"),
  sep="\t",
  row.names = FALSE,
  quote = FALSE
)


### pcr only: redo with "allreads" later. empty
writeMageckOut(
  "14d_pcr",
  list_14d_pcr_normoxia,
  list_lib_pcr
)

writeMageckOut(
  "3d_pcr",
  list_3d_pcr_normoxia,
  list_lib_pcr
)




writeMageckOut(
  "14d_mip_dedup",
  list_14d_mip_dedup_normoxia,
  list_lib_mip_dedup
)

writeMageckOut( 
  "3d_mip_dedup",
  list_3d_mip_dedup_normoxia,
  list_lib_mip_dedup
)


writeMageckOut( 
  "14d_mip_nodedup",
  list_14d_mip_nodedup_normoxia,
  list_lib_mip_nodedup
)

writeMageckOut( 
  "3d_mip_nodedup",
  list_3d_mip_nodedup_normoxia,
  list_lib_mip_nodedup
)






################################################################################
######### Efficiency of MIP capture ############################################
################################################################################


dir_subsamp <- "/corgi/otherdataset/crispr_padlock/subsamp"

out_summary <- NULL
#listsum <- list.files(dir_subsamp, pattern = "*15000000#0")
listsum <- list.files(dir_subsamp, pattern = "*allreads#0")
for(f in listsum){
  curcond <- str_split_fixed(f,"#",3)[1]
  numread <- as.integer(str_split_fixed(f,"#",3)[2])  #if NA, then this is max reads
  print(f)
  
  dat <- read.csv(file.path(dir_subsamp, f),sep=",")
  
  df <- data.frame(
    cond=curcond,
    num_umi=sum(dat$cnt)
  )
  out_summary <- rbind(df, out_summary)
}
#out_summary[str_detect(out_summary$cond,"P29968_1001"),] 


out_summary$filename <- paste0(out_summary$cond,".csv")

#out_summary$filename <- out_summary$cond 
#out_summary$filename[str_detect(out_summary$filename,"dedup")] <- paste0(out_summary$filename[str_detect(out_summary$filename,"dedup")],".csv")
#out_summary$filename[str_detect(out_summary$filename,"nodedup")] <- paste0(out_summary$filename[str_detect(out_summary$filename,"nodedup")],".csv")
#out_summary$filename[!str_detect(out_summary$filename,"dedup")] <- paste0(out_summary$filename[!str_detect(out_summary$filename,"dedup") & !str_detect(out_summary$filename,"dedup")],"_dedup.csv")

head(samplemeta_padlock_nodedup)

####################### Analyze for 3d and 14d sample
out_summary_sub <- merge(samplemeta_padlock_dedup, out_summary)

out_summary_sub <- out_summary_sub[,c("samplename","time","cells","protocol","dedup","num_umi")]
out_summary_sub$eff <- out_summary_sub$num_umi/out_summary_sub$cells

out_summary_sub[out_summary_sub$time %in% c("3d","14d"),]
#13%  ... with all reads!
#3.8% with 100k reads
#8.6% with 1M reads
#12.1% with 10M reads
#13.0% with 15M reads
#14.6% with all reads reads

mean(out_summary[out_summary$time %in% c("3d","14d"),]$eff)




####################### Analyze for digested DNA sample
out_summary_sub <- merge(samplemeta, out_summary)

out_summary_sub <- out_summary_sub[,c("filename","samplename","time","cells","protocol","dedup","num_umi","treatment")]
out_summary_sub$cells <- 150e3 ## error in samplemeta
out_summary_sub$eff <- out_summary_sub$num_umi/out_summary_sub$cells

out_summary_sub[out_summary_sub$treatment %in% c("digested_gDNA","non-digested_gDNA") & out_summary_sub$dedup,]
#13%  ... with all reads!
#3.8% with 100k reads
#8.6% with 1M reads
#12.1% with 10M reads
#13.0% with 15M reads
#14.6% with all reads reads

mean(out_summary[out_summary$time %in% c("3d","14d"),]$eff)
