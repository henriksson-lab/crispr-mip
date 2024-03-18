
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

  writeLines(all_cmd_mageck_mle, file.path(dir_mageck, paste0(comparison_name, ".mle.sh")))
  writeLines(all_cmd_mageck, file.path(dir_mageck, paste0(comparison_name, ".mageck.sh")))
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



