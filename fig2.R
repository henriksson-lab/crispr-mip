
################################################################################
################## Figure 2c - Theoretical UMI distribution ####################
################################################################################

# python code to get one sgRNA: umi_dedup_stats.ipynb

#1% error rate from literature; 3.5 fitted here

errorrate <- 0.01
umilen <- 13

if(FALSE){
  df <- data.frame(
    num_mismatch=0:5
  )
  df$p <- dpois(df$num_mismatch,umilen*errorrate)
  ggplot(df, aes(num_mismatch, p)) + geom_bar(stat = "identity")

}





umi_dist <- read.csv("umi_dist.csv") ######## From python
cnt_main_umi <- umi_dist$cnt[1]
cnt_total_reads <- sum(umi_dist$cnt)

#Compute empirical distribution
#How many of the reads are this umi?
dist_cnt <- sqldf::sqldf("select sum(cnt) as cnt, dist from umi_dist group by dist")
dist_cnt_empirical <- dist_cnt
dist_cnt <- merge(dist_cnt, data.frame(dist=0:13), all=TRUE)
dist_cnt[is.na(dist_cnt)] <- 0
dist_cnt$cnt <- dist_cnt$cnt/sum(dist_cnt$cnt)
dist_cnt$src <- "Empirical"
estimate_frac_thisumi <- sum(dist_cnt$cnt[dist_cnt$dist<6]) #/ sum(dist_cnt$cnt) ######## This number used later

ggplot(dist_cnt, aes(dist, cnt, fill=src)) + geom_bar(stat = "identity") #For testing


### Estimate theoretical distance to one other UMI. Because of central limit theorem, will we expect a more centralized distribution?
dist_bin_other <- data.frame(
  dist=0:13, 
  cnt=dbinom(0:13, prob=0.75, size=13)*(1-estimate_frac_thisumi),
  src="Other UMIs"
)



#Fit lambda
x.pois <- c()
for(i in 1:6){
  print(i)
  x.pois <- c(x.pois, rep(dist_cnt_empirical$dist[i], times=dist_cnt_empirical$cnt[i]))
}
my.mle <- MASS::fitdistr(x.pois, densfun="poisson")
fit_lambda <- my.mle$estimate
fit_lambda ###### what does this correspond to in terms of sequencing errors? Actually, should be a poisson distribution too I think?


### Estimate theoretical distance to the same UMI, but with sequencing errors
dist_bin_pois <- data.frame(
  dist=0:13, 
  cnt=dpois(0:13, lambda = fit_lambda)*estimate_frac_thisumi,
  src="This UMI (poisson)"
)
ggplot(rbind(dist_cnt,dist_bin_pois,dist_bin_other), aes(dist, cnt, fill=src)) + geom_bar(stat = "identity", position = "dodge2")  #For testing



###### binomial: mean is n*p; n=13
dist_cnt_this <- dist_cnt_empirical
dist_cnt_this$cnt <- dist_cnt_this$cnt / sum(dist_cnt_this$cnt)
mean_this <- sum(dist_cnt_this$cnt[dist_cnt_this$dist<6] * dist_cnt_this$dist[dist_cnt_this$dist<6]) / sum(dist_cnt_this$cnt)
est_this_prob <- mean_this/13
est_this_prob <- 0.035

### Estimate theoretical distance to the same UMI, but with sequencing errors ----------- new
dist_bin_this <- data.frame(
    dist=0:13, 
    cnt=dbinom(0:13, prob=est_this_prob, size=13)*estimate_frac_thisumi,
    src="This UMI"
)
ggplot(rbind(
  dist_cnt,
  dist_bin_this,
#  dist_bin_pois,
  dist_bin_other
), aes(dist, cnt, fill=src)) + 
  geom_bar(stat = "identity", position = "dodge2") +
  theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(name = "Difference (bp)", breaks=c(0:13))+#) +
  ylab("Frequency")
ggsave("out/fig2c.svg", width = 5, height = 2)


#ggplot(dist_cnt, aes(dist, cnt, fill=src)) + geom_bar(stat = "identity", position = "dodge2")






################################################################################
############ Figure 2d  / sup fig 1 - How to best deduplicate ##################
################################################################################


# Generated from: simulate_UMI.ipynb
dat <- read.csv("out/umitools_stat.csv")
dat <- sqldf::sqldf("select avg(detected_umi) as detected_umi, num_umi, mincount, threshold, num_read_per_umi from dat group by num_umi, mincount, threshold, num_read_per_umi")

dat$error <- (dat$detected_umi - dat$num_umi)/dat$num_umi
range(dat$error)

list_plot <- list()
for(num_read_per_umi in unique(dat$num_read_per_umi)){
  for(mincount in unique(dat$mincount)){
    breaks <- c(-1, 0, 1,2,3,4) #cannot get below -1
    mycolours <- c("red", "white", "blue")
    yvals <- sort(unique(dat$threshold))
    
    list_plot[[paste("mc",mincount, "rpu",num_read_per_umi)]] <-
      ggplot(dat[dat$mincount==mincount,], aes(num_umi, threshold, fill=error)) + 
      scale_x_log10() +
      geom_tile() + 
      #    scale_fill_gradientn(colors = mycolours, breaks=breaks,labels=format(breaks))+
      scale_fill_gradient2() + 
      
      #scale_fill_gradientn(limits = c(0.1,100000), breaks=c(0.1,1,10,100,1000,10000,100000), labels=c(0.1,1,10,100,1000,10000,100000),
       #                    colours=c("navyblue", "gray", "darkorange1"))+
      
      xlab(paste("#UMI, thr =",mincount,"Read/UMI",num_read_per_umi))+
      scale_y_continuous(breaks=yvals) + 
      theme_gray(base_size = 14) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  } 
}
ptot <- do.call(egg::ggarrange,list_plot)
ptot
ggsave("out/supfig1.svg", ptot, width = 20, height = 10)




################################################################################
################## Figure UNUSED - sgrna coverage analysis #####################
################################################################################

if(FALSE){

  library(stringr)
  library(dplyr)
  
  # Loop over all files, collect depth, 
  grna_cnt <- data.frame(
    f=list.files(countdir),
    num_grna=NA
  )
  for(i in seq_along(grna_cnt$f)){
    thefile <- file.path(countdir,grna_cnt$f[i])
    onecnt <- read.csv(thefile)
    grna_cnt$num_grna[i] <- sum(onecnt$cnt>0)
  }
  grna_cnt$rep <- as.integer(str_split_fixed(grna_cnt$f,"#",3)[,3])
  grna_cnt$numread <- as.integer(str_split_fixed(grna_cnt$f,"#",3)[,2])
  grna_cnt$filename <- paste0(str_split_fixed(grna_cnt$f,"#",3)[,1],".csv")
  
  
  ################
  samplemeta$group_cnt <- NA
  
  samplemeta$group_cnt[samplemeta$filename %in% list_3d_pcr_normoxia] <- "list_3d_pcr_normoxia"
  samplemeta$group_cnt[samplemeta$filename %in% list_3d_mip_dedup_normoxia] <- "list_3d_mip_dedup_normoxia"
  samplemeta$group_cnt[samplemeta$filename %in% list_3d_mip_nodedup_normoxia] <- "list_3d_mip_nodedup_normoxia"
  
  samplemeta$group_cnt[samplemeta$filename %in% list_14d_pcr_normoxia] <- "list_14d_pcr_normoxia"
  samplemeta$group_cnt[samplemeta$filename %in% list_14d_mip_dedup_normoxia] <- "list_14d_mip_dedup_normoxia"
  samplemeta$group_cnt[samplemeta$filename %in% list_14d_mip_nodedup_normoxia] <- "list_14d_mip_nodedup_normoxia"
  
  samplemeta$group_cnt[samplemeta$filename %in% list_lib_pcr] <- "list_lib_pcr"
  samplemeta$group_cnt[samplemeta$filename %in% list_lib_mip_dedup] <- "list_lib_mip_dedup"
  samplemeta$group_cnt[samplemeta$filename %in% list_lib_mip_nodedup] <- "list_lib_mip_nodedup"
  
  grna_cnt <- merge(grna_cnt, samplemeta[,c("filename","group_cnt")])
  
  
  #######
  # simple plot of each point
  plot_grna_count_with_confidence <- function(grna_sub){
    ggplot(grna_sub, aes(log10(numread), num_grna, color=group_cnt)) + geom_point() + xlab("#reads")
  }
  
  
  ##########
  # plot with conf interval
  plot_grna_count_with_confidence <- function(grna_sub){
    grouped <- group_by(grna_sub, group_cnt, numread)
    forgg <- summarise(grouped, mean=mean(num_grna), sd=sd(num_grna))
    ggplot(data = forgg, aes(x = log10(numread), group = group_cnt)) + 
      geom_line(aes(y = mean, color = group_cnt), size = 1) + 
      geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = group_cnt), alpha = .2) +
      xlab("#Reads") + ylab("#sgRNA")
  }
  
  
  grna_sub <- grna_cnt[!is.na(grna_cnt$group_cnt) & str_detect(grna_cnt$group_cnt,"lib"),]
  p1 <- plot_grna_count_with_confidence(grna_sub) + xlab("#Reads, plasmid")
  
  grna_sub <- grna_cnt[!is.na(grna_cnt$group_cnt) & str_detect(grna_cnt$group_cnt,"14d") & str_detect(grna_cnt$group_cnt,"normoxia"),]
  p2 <- plot_grna_count_with_confidence(grna_sub) + xlab("#Reads, genomic DNA (14d dropout)")
  
  egg::ggarrange(p1,p2,ncol=2)
  
  
}





################################################################################
############### Compute what max reads is for each lib #########################
################################################################################



maxdepth_per_cond <- NULL
listsum <- list.files(dir_mageck, pattern = "*count")
for(f in listsum){
  if(str_detect(f,"allreads")){
    curcond <- str_split_fixed(f,"#",3)[1]
    dat <- read.csv(file.path(dir_mageck, f),sep="\t")
    df <- data.frame(
      cond=curcond,
      maxreads=exp(mean(log(
        colSums(dat[,-c(1:2)])
      )))
    )
    maxdepth_per_cond <- rbind(df, maxdepth_per_cond) 
  }
}

################################################################################
################## Figure 1efg - screen hits vs depth ##########################
################################################################################

#must have run for_mageck.R, and mageck, before this section

dir_mageck <- "/corgi/otherdataset/crispr_padlock/for_mageck"


################### plotting: mageck TEST analysis  - sgrna

min_read_cutoff <- 10000

out_summary <- NULL
listsum <- list.files(dir_mageck, pattern = "*sgrna_summary.txt")
for(f in listsum){
  if(!str_detect(f,"mageckmle")){
    curcond <- str_split_fixed(f,"#",3)[1]
    numread <- as.integer(str_split_fixed(f,"#",3)[2])  #if NA, then this is max reads
    
    dat <- read.csv(file.path(dir_mageck, f),sep="\t")
    numsig <- sum(dat$FDR<1e-5)
    
    df <- data.frame(
      cond=curcond,
      numread=numread,
      numsig=numsig
    )
    out_summary <- rbind(df, out_summary) 
  }
}
out_summary
#out_summary$numread[is.na(out_summary$numread)] <- 1e8  ## set for max reads



if(FALSE){
  #Not included, because hard to set a good x-value for max reads
  out_summary_maxread <- out_summary[is.na(out_summary$numread),]
  out_summary_maxread <- merge(maxdepth_per_cond,out_summary_maxread)
  out_summary_maxread$numread <- max(out_summary$numread, na.rm = TRUE)*1.5 #out_summary_maxread$maxreads
  out_summary_maxread$sd <- 0
  out_summary_maxread$mean <- out_summary_maxread$numread
}

########## plot with conf interval - d14
grouped <- group_by(out_summary[!is.na(out_summary$numread) & out_summary$numread>min_read_cutoff,], cond, numread)
forgg <- summarise(grouped, mean=mean(numsig), sd=sd(numsig))
p_grna_test_14d <- ggplot(data = forgg[str_detect(forgg$cond,"14"),], aes(x = numread, group = cond)) + 
  geom_line(aes(y = mean, color = cond), size = 1) + 
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = cond), alpha = .2) +
  xlab("#reads") + ylab("MAGeCK TEST #sgRNA (d14)") + 
  scale_x_log10() #+
  #geom_point(data = out_summary_maxread[str_starts(out_summary_maxread$cond,"14d"),], aes(x=numread, y=numsig, color = cond))

p_grna_test_14d



########## plot with conf interval - d3
grouped <- group_by(out_summary[!is.na(out_summary$numread) & out_summary$numread>min_read_cutoff,], cond, numread)
forgg <- summarise(grouped, mean=mean(numsig), sd=sd(numsig))
p_grna_test_3d <- ggplot(data = forgg[str_detect(forgg$cond,"3"),], aes(x = numread, group = cond)) + 
  geom_line(aes(y = mean, color = cond), size = 1) + 
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = cond), alpha = .2) +
  xlab("#reads") + ylab("MAGeCK TEST #sgRNA (d3)") + scale_x_log10()





################### plotting: mageck TEST analysis  - genes

out_summary <- NULL
listsum <- list.files(dir_mageck, pattern = "*gene_summary.txt")
for(f in listsum){
  if(!str_detect(f,"mageckmle")){
    curcond <- str_split_fixed(f,"#",3)[1]
    numread <- as.integer(str_split_fixed(f,"#",3)[2])
    
    dat <- read.csv(file.path(dir_mageck, f),sep="\t")
    numsig <- sum(dat$neg.fdr<0.05)

    df <- data.frame(
      cond=curcond,
      numread=numread,
      numsig=numsig
    )
    out_summary <- rbind(df, out_summary)
  }
}
#out_summary$numread[is.na(out_summary$numread)] <- 1e8  ## set for max reads


########## plot with conf interval - d3
grouped <- group_by(out_summary[!is.na(out_summary$numread) & out_summary$numread>min_read_cutoff,], cond, numread)
forgg <- summarise(grouped, mean=mean(numsig), sd=sd(numsig))
p_gene_test_3d <- ggplot(data = forgg[str_detect(forgg$cond,"3"),], aes(x = numread, group = cond)) + 
  geom_line(aes(y = mean, color = cond), size = 1) + 
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = cond), alpha = .2) +
  xlab("#reads") + ylab("MAGeCK TEST #gene (d3)") + scale_x_log10()

########## plot with conf interval - d14
p_gene_test_14d <- ggplot(data = forgg[str_detect(forgg$cond,"14"),], aes(x = numread, group = cond)) + 
  geom_line(aes(y = mean, color = cond), size = 1) + 
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = cond), alpha = .2) +
  xlab("#reads") + ylab("MAGeCK TEST #gene (d14)") + scale_x_log10()

p_gene_test_14d



################### plotting: mageck MLE analysis - genes

out_summary <- NULL
listsum <- list.files(dir_mageck, pattern = "*.mageckmle.gene_summary.txt")
for(f in listsum){
  curcond <- str_split_fixed(f,"#",3)[1]
  numread <- as.integer(str_split_fixed(f,"#",3)[2])
  
  dat <- read.csv(file.path(dir_mageck, f),sep="\t")
  numsig <- sum(dat$eff.fdr<0.05)
  
  df <- data.frame(
    cond=curcond,
    numread=numread,
    numsig=numsig
  )
  out_summary <- rbind(df, out_summary)
}
#out_summary$numread[is.na(out_summary$numread)] <- 1e8  ## set for max reads


########## plot with conf interval
grouped <- group_by(out_summary[!is.na(out_summary$numread) & out_summary$numread>min_read_cutoff,], cond, numread)
forgg <- summarise(grouped, mean=mean(numsig), sd=sd(numsig))
p_gene_mle_3d <- ggplot(data = forgg[str_detect(forgg$cond,"3"),], aes(x = numread, group = cond)) + 
  geom_line(aes(y = mean, color = cond), size = 1) + 
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = cond), alpha = .2) +
  xlab("#reads") + ylab("MAGeCK MLE #gene (d3)") + scale_x_log10()

p_gene_mle_14d <- ggplot(data = forgg[str_detect(forgg$cond,"14"),], aes(x = numread, group = cond)) + 
  geom_line(aes(y = mean, color = cond), size = 1) + 
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = cond), alpha = .2) +
  xlab("#reads") + ylab("MAGeCK MLE #gene (d14)") + scale_x_log10()

########################### Figure 1efg and supplemental fig2 --  assemble it all

format_one_panel <- function(p){
  p + theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    scale_x_log10(limits = c(1e4,2e7)) 
}

### Supplemental figure
ptot <- egg::ggarrange(
  format_one_panel(p_gene_mle_3d)+ylim(0,200), format_one_panel(p_gene_mle_14d)+ylim(0,200),
  format_one_panel(p_gene_test_3d)+ylim(0,100), format_one_panel(p_gene_test_14d)+ylim(0,100),
  format_one_panel(p_grna_test_3d)+ylim(0,1200), format_one_panel(p_grna_test_14d)+ylim(0,1200)
)
ptot
ggsave(filename = "out/supfig3.svg", ptot, width = 6*2, height = 3*3)

### For fig1efg
ptot <- egg::ggarrange(
  format_one_panel(p_grna_test_14d)+ylim(0,2100),
  format_one_panel(p_gene_test_14d)+ylim(0,130),
  format_one_panel(p_gene_mle_14d)+ylim(0,130),
  ncol = 3
)
ptot
ggsave(filename = "out/fig1_efg.svg", ptot, width = 6*2*1.5, height = 2*1.5)




################################################################################
################## Export as XLS ###############################################
################################################################################



if(FALSE){
  install.packages("xlsx")
}

library("xlsx")

dir_mageck <- "/corgi/otherdataset/crispr_padlock/for_mageck"


list_cond <- c(
  "PCR MAGeCK test",
  "PCR MAGeCK MLE",
  "CRISPR-MIP dedup MAGeCK test",
  "CRISPR-MIP dedup MAGeCK mle"
)
list_fname <- c(
  "14d_pcr#allreads#0.gene_summary.txt",
  "14d_pcr#allreads#0.mageckmle.gene_summary.txt",
  "14d_mip_dedup#allreads#0.gene_summary.txt",
  "14d_mip_dedup#allreads#0.mageckmle.gene_summary.txt"
)
for(i in 1:length(list_cond)){
  dat <- read.csv(file.path(dir_mageck, list_fname[i]),sep="\t", as.is = TRUE, check.names = FALSE)
  if("neg|fdr" %in% colnames(dat)){
    numsig <- sum(dat$`neg|fdr`<0.05)
  } else {
    numsig <- sum(dat$`eff|fdr`<0.05)
    dat <- dat[order(dat$`eff|fdr`),]
    rownames(dat) <- NULL
  }
  print(paste(list_cond[i], "---", list_fname[i], "--- num sig:", numsig))
  write.xlsx(dat, file = "/corgi/otherdataset/crispr_padlock/out_sup_table1.xlsx",
             sheetName = list_cond[i], append = i!=1, row.names = FALSE)  
}





################################################################################
######### Behavior of mageck for PCR -- MLE ####################################
################################################################################


dir_mageck2 <- "/corgi/otherdataset/crispr_padlock/for_mageck"

#min_read_cutoff <- 10000

out_summary <- list()
listsum <- list.files(dir_mageck2, pattern = "*.mageckmle.gene_summary.txt")
for(f in listsum){
  #print(f)
  if(str_detect(f,fixed("#0.mageckmle.gene_summary.txt")) & str_detect(f,fixed("14d_"))){
    print(f)
    curcond <- str_split_fixed(f,"#",3)[1]
    #numread <- as.integer(str_split_fixed(f,"#",3)[2])  #if NA, then this is max reads
    
    dat <- read.csv(file.path(dir_mageck2, f),sep="\t")
    dat$depth <- as.integer(str_split_i(f,fixed("#"),2))
    dat$cond <- curcond
    
    out_summary[[f]] <- dat
  }
}
out_summary <- do.call(rbind, out_summary)
out_summary <- out_summary[out_summary$depth>=1e4,]

#mat_mle_effbeta <- acast(out_summary[,c("Gene","depth","eff.beta")],Gene~depth, value.var="eff.beta")
p1 <- ggplot(out_summary[str_detect(out_summary$cond,"14d_pcr"),], aes(depth, eff.beta, color=Gene)) + 
  geom_line() + 
  scale_x_log10() + 
  theme(legend.position = "none") +
  ylab("eff.beta, PCR")
#  theme_bw()
#  xlab()
p2 <- ggplot(out_summary[str_detect(out_summary$cond,"14d_mip_dedup"),], aes(depth, eff.beta, color=Gene)) + 
  geom_line() + 
  scale_x_log10() + 
  theme(legend.position = "none") +
  ylab("eff.beta, MIP dedup")
p3 <- ggplot(out_summary[str_detect(out_summary$cond,"14d_mip_nodedup"),], aes(depth, eff.beta, color=Gene)) + 
  geom_line() + 
  scale_x_log10() + 
  theme(legend.position = "none") +
  ylab("eff.beta, MIP no dedup")
egg::ggarrange(p1,p2,p3)


if(FALSE){
  ggplot(out_summary[str_detect(out_summary$cond,"14d_pcr"),], aes(depth, (eff.p.value), color=Gene)) + 
    geom_line() + 
    theme(legend.position = "none")
  
  ggplot(out_summary[str_detect(out_summary$cond,"14d_pcr"),], aes(depth, (eff.fdr), color=Gene)) + 
    geom_line() + 
    theme(legend.position = "none")
  
  ggplot(out_summary[str_detect(out_summary$cond,"14d_pcr"),], aes(depth, (eff.fdr), color=Gene)) + 
    geom_line() + 
    theme(legend.position = "none")
  
}




################################################################################
######### Behavior of mageck for PCR -- test ###################################
################################################################################


dir_mageck2 <- "/corgi/otherdataset/crispr_padlock/for_mageck"

out_summary <- list()
listsum <- list.files(dir_mageck2, pattern = "*#0.sgrna_summary.txt")
for(f in listsum){
  #print(f)
  if(!str_detect(f,fixed("mageckmle")) & str_detect(f,fixed("14d_"))){
    print(f)
    curcond <- str_split_fixed(f,"#",3)[1]
    #numread <- as.integer(str_split_fixed(f,"#",3)[2])  #if NA, then this is max reads
    
    dat <- read.csv(file.path(dir_mageck2, f),sep="\t")
    dat$depth <- as.integer(str_split_i(f,fixed("#"),2))
    dat$cond <- curcond
    
    out_summary[[f]] <- dat
  }
}
out_summary <- do.call(rbind, out_summary)


if(FALSE){
  p1 <- ggplot(out_summary[str_detect(out_summary$cond,"14d_pcr"),], aes(depth, FDR, color=sgrna)) + 
    geom_line() + 
    scale_x_log10() + 
    theme(legend.position = "none") +
    ylab("FDR, PCR")
  p2 <- ggplot(out_summary[str_detect(out_summary$cond,"14d_mip_dedup"),], aes(depth, FDR, color=sgrna)) + 
    geom_line() + 
    scale_x_log10() + 
    theme(legend.position = "none") +
    ylab("FDR, MIP dedup")
  p3 <- ggplot(out_summary[str_detect(out_summary$cond,"14d_mip_nodedup"),], aes(depth, FDR, color=sgrna)) + 
    geom_line() + 
    scale_x_log10() + 
    theme(legend.position = "none") +
    ylab("FDR, MIP no dedup")
  egg::ggarrange(p1,p2,p3)
}







#mat_mle_effbeta <- acast(out_summary[,c("Gene","depth","eff.beta")],Gene~depth, value.var="eff.beta")
##FDR scores
if(FALSE){
  ggplot(out_summary[str_detect(out_summary$cond,"14d_pcr"),], aes(depth, FDR, color=sgrna)) + 
    geom_line() + 
    scale_x_log10()+
    theme(legend.position = "none")
  
}

ggplot(out_summary[str_detect(out_summary$cond,"14d_pcr"),], aes(depth, FDR)) + 
  geom_density_2d_filled(alpha = 0.5)+
  #geom_line() + 
  scale_x_log10()+
  theme(legend.position = "none")


######### per sgRNA lfc it stabilizes
ggplot(out_summary[str_detect(out_summary$cond,"14d_pcr"),], aes(depth, LFC, color=sgrna)) + 
  geom_line() + 
  scale_x_log10()+
  theme(legend.position = "none")


### Variance of control, normalized by mean ---  appears somewhat stable, but does increase(?)
ggplot(out_summary[str_detect(out_summary$cond,"14d_pcr"),], aes(depth, control_var/control_mean, color=sgrna)) + 
  geom_line() + 
  scale_x_log10()+
  scale_y_log10()+
  theme(legend.position = "none")
ggplot(out_summary[str_detect(out_summary$cond,"14d_pcr"),], aes(depth, control_var/control_mean)) + 
  geom_density_2d_filled(alpha = 0.5)+
  scale_x_log10()+
  scale_y_log10()+
  theme(legend.position = "none")


### Adjusted variance of control/mean, goes to hell.
#eq 4 in https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0554-4   is the adjusted_mean from NB assumption
ggplot(out_summary[str_detect(out_summary$cond,"14d_pcr"),], aes(depth, adj_var/control_mean, color=sgrna)) + 
  geom_line() + 
  scale_x_log10()+
  scale_y_log10()+
  theme(legend.position = "none")


ggplot(out_summary[str_detect(out_summary$cond,"14d_mip_dedup"),], aes(depth, adj_var/control_mean, color=sgrna)) + 
  geom_line() + 
  scale_x_log10()+
  scale_y_log10()+
  theme(legend.position = "none")

ggplot(out_summary[str_detect(out_summary$cond,"14d_mip_nodedup"),], aes(depth, adj_var/control_mean, color=sgrna)) + 
  geom_line() + 
  scale_x_log10()+
  scale_y_log10()+
  theme(legend.position = "none")

# 
# 
# out_summary[out_summary$score>250,]
# 
# out_summary$ourz <-  log2(out_summary$treat_mean/out_summary$control_mean)
# ggplot(out_summary, aes(depth, ourz, color=sgrna)) + 
#   geom_line() + 
#   scale_x_log10()+
#   theme(legend.position = "none")
# 
# ggplot(out_summary[out_summary$Gene %in% c("CDK1"),], aes(depth, ourz, color=sgrna)) + 
#   geom_line() + 
#   scale_x_log10()+
#   theme(legend.position = "none")


#############
############# per gene analysis -- but the issue is on the level of sgRNAs
#############


out_summary <- list()
listsum <- list.files(dir_mageck2, pattern = "*#0.gene_summary.txt")
for(f in listsum){
  #print(f)
  if(!str_detect(f,fixed("mageckmle")) & str_detect(f,fixed("14d_pcr"))){
    print(f)
    curcond <- str_split_fixed(f,"#",3)[1]
    #numread <- as.integer(str_split_fixed(f,"#",3)[2])  #if NA, then this is max reads
    
    dat <- read.csv(file.path(dir_mageck2, f),sep="\t")
    dat$depth <- as.integer(str_split_i(f,fixed("#"),2))
    
    out_summary[[f]] <- dat
  }
}
out_summary <- do.call(rbind, out_summary)
out_summary

#####  this looks fine
ggplot(out_summary, aes(depth, pos.lfc, color=id)) + 
  geom_line() + 
  scale_x_log10()+
  theme(legend.position = "none")

ggplot(out_summary, aes(depth, pos.lfc, color=id)) + 
  geom_line() + 
  scale_x_log10()+
  theme(legend.position = "none")




################################################################################
######### Test NB assumption over PCR data #####################################
################################################################################




#mageck test -k  /corgi/otherdataset/crispr_padlock/for_mageck/14d_pcr#allreads#0.count -t 0,1,2 -c 3,4,5,6 --control-sgrna /corgi/otherdataset/crispr_padlock/for_mageck/14d_pcr#allreads#0.control -n /corgi/otherdataset/crispr_padlock/for_mageck/14d_pcr#allreads#0
#mageck test -k  /corgi/otherdataset/crispr_padlock/for_mageck/14d_mip_dedup#allreads#0.count -t 0,1,2 -c 3,4,5 --control-sgrna /corgi/otherdataset/crispr_padlock/for_mageck/14d_mip_dedup#allreads#0.control -n /corgi/otherdataset/crispr_padlock/for_mageck/14d_mip_dedup#allreads#0


#mageck test -k  /corgi/otherdataset/crispr_padlock/for_mageck/14d_mip_dedup#10000000#0.count -t 0,1,2 -c 3,4,5 --control-sgrna /corgi/otherdataset/crispr_padlock/for_mageck/14d_mip_dedup#10000000#0.control -n /corgi/otherdataset/crispr_padlock/for_mageck/14d_mip_dedup#10000000#0
#mageck test -k  /corgi/otherdataset/crispr_padlock/for_mageck/14d_mip_dedup#100000#0.count -t 0,1,2 -c 3,4,5 --control-sgrna /corgi/otherdataset/crispr_padlock/for_mageck/14d_mip_dedup#100000#0.control -n /corgi/otherdataset/crispr_padlock/for_mageck/14d_mip_dedup#100000#0


#mageck test -k  /corgi/otherdataset/crispr_padlock/for_mageck/14d_mip_nodedup#10000000#0.count -t 0,1,2 -c 3,4,5 --control-sgrna /corgi/otherdataset/crispr_padlock/for_mageck/14d_mip_nodedup#10000000#0.control -n /corgi/otherdataset/crispr_padlock/for_mageck/14d_mip_nodedup#10000000#0
#mageck test -k  /corgi/otherdataset/crispr_padlock/for_mageck/14d_mip_nodedup#100000#0.count -t 0,1,2 -c 3,4,5 --control-sgrna /corgi/otherdataset/crispr_padlock/for_mageck/14d_mip_nodedup#100000#0.control -n /corgi/otherdataset/crispr_padlock/for_mageck/14d_mip_nodedup#100000#0

gcproc <- function(grna){
  (stringr::str_count(grna, "G") + stringr::str_count(grna, "C")) / str_length(grna)
}

readx <- function(f) as.double(str_split(str_remove_all(str_remove_all(readLines(f),fixed("[")),fixed("]")),", ")[[1]])
readsgrna <- function(f) str_split(str_remove_all(str_remove_all(str_remove_all(readLines(f),fixed("'")),fixed("[")),fixed("]")),", ")[[1]]

load_meanvar <- function(dir){
  data.frame(
    sgrna=readsgrna(file.path(dir,"all_sgrna.txt")),
    x=readx(file.path(dir,"all_x.txt")),
    y=readx(file.path(dir,"all_y.txt")),
    w=readx(file.path(dir,"all_w.txt"))
  )
}

df_mip_dedup_10m  <- load_meanvar("/corgi/otherdataset/crispr_padlock/meanvar_relation/dedup/10000000")
df_mip_dedup_100k <- load_meanvar("/corgi/otherdataset/crispr_padlock/meanvar_relation/dedup/100000")
dim(df_mip_dedup_10m)
dim(df_mip_dedup_100k)

df_mip_nodedup_10m  <- load_meanvar("/corgi/otherdataset/crispr_padlock/meanvar_relation/nodedup/10000000")
df_mip_nodedup_100k <- load_meanvar("/corgi/otherdataset/crispr_padlock/meanvar_relation/nodedup/100000")
dim(df_mip_nodedup_10m)
dim(df_mip_nodedup_100k)

if(FALSE){
  #Thr=1 makes no difference
  df_mip_dedup_10m  <- load_meanvar("/corgi/otherdataset/crispr_padlock/meanvar_relation/thr1_dedup/10000000")
  df_mip_dedup_100k <- load_meanvar("/corgi/otherdataset/crispr_padlock/meanvar_relation/thr1_dedup/100000")
  df_mip_nodedup_10m  <- load_meanvar("/corgi/otherdataset/crispr_padlock/meanvar_relation/thr1_nodedup/10000000")
  df_mip_nodedup_100k <- load_meanvar("/corgi/otherdataset/crispr_padlock/meanvar_relation/thr1_nodedup/100000")
}


df_pcr_10m  <- load_meanvar("/corgi/otherdataset/crispr_padlock/meanvar_relation/pcr/10000000")
df_pcr_100k <- load_meanvar("/corgi/otherdataset/crispr_padlock/meanvar_relation/pcr/100000")
dim(df_pcr_10m)
dim(df_pcr_100k)




make_meanvar_plot <- function(df){
  df$gc <- gcproc(df$sgrna)
  
  #df <- df[df$gc<0.6 & df$gc>0.4,] #does not improve fit
  
  themod <- lm(y~x, data = df, weights = df$w)  #consistent with mageck calculations
  ggplot(df, aes(x,y)) + 
    geom_point() + 
    geom_abline(slope=themod$coefficients[2], intercept = themod$coefficients[1], color="red")
  
  #if k<1:  k=1
  #if b<0:  b=0   #this is mageck adjustments
  ggplot(df, aes(x,y,color=gc)) + 
    geom_point() + 
    geom_abline(slope=themod$coefficients[2], intercept = themod$coefficients[1], color="red") +
    geom_abline(slope=max(themod$coefficients[2],1), intercept = max(themod$coefficients[1],0), color="blue") +
    geom_abline(slope=1, intercept = 0, color="gray") +
    xlab("log2(1+mean)")+
    ylab("log2(1+var)")+
    theme_bw()+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    xlim(0,15)+
    ylim(0,17)
}

ptot <- egg::ggarrange(
  make_meanvar_plot(df_pcr_100k)+ggtitle("pcr 100k"),
  make_meanvar_plot(df_mip_nodedup_100k)+ggtitle("mip no dedup 100k"),
  make_meanvar_plot(df_mip_dedup_100k)+ggtitle("mip dedup 100k"),
  
  make_meanvar_plot(df_pcr_10m)+ggtitle("pcr 10m"),
  make_meanvar_plot(df_mip_nodedup_10m)+ggtitle("mip no dedup 10m"),
  make_meanvar_plot(df_mip_dedup_10m)+ggtitle("mip dedup 10m"),
  nrow=2
)
ptot
ggsave(plot=ptot, "newout/mageck_test_meanvar.svg", width = 8, height = 4)


#TODO: is dedup for 10M too aggressive???

###### Can we fit an NB distribution?
#logVar = logMean + log(1+mean/kappa)



################################################################################
####### Supplemental figure xxx: Simulate mean-var of poisson ##################
################################################################################



make_plot_sim <- function(for_lambda, npoint, nsamp) {
  dat <- matrix(rpois(n = nsamp*npoint, lambda = for_lambda), nrow = nsamp, ncol = npoint)
  
  df <- data.frame(
    m=colMeans(dat),
    v=colVars(dat, useNames = FALSE)  #note: this uses n/(n−1)∗mean((x−center)**2)  ; same as getVars in mageck, mageckMathFunc.py
  )
  ggplot(df, aes(log2(1+m), log2(1+v))) +
    geom_point() +
    xlim(0,log2(for_lambda*2)) +
    ylim(0,log2(for_lambda*2)) +
    geom_abline(slope=1, color="gray") + 
    theme_bw()+
    xlab("log2(1+mean)")+
    ylab("log2(1+var)")+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())  
}


p1 <- make_plot_sim(
  for_lambda = 2**5,
  npoint = 500,
  nsamp = 6
)+xlim(0,10) + ylim(0,10) 
p2 <- make_plot_sim(
  for_lambda = 2**10,
  npoint = 500,
  nsamp = 6
)+xlim(0,10) + ylim(0,10) 


ptot <- egg::ggarrange(p1,p2, nrow=1)
ggsave(plot = ptot, "sim fitted poiss_n6.pdf", width = 4, height = 2)




p1 <- make_plot_sim(
  for_lambda = 2**5,
  npoint = 500,
  nsamp = 3
)+xlim(0,10) + ylim(0,10) 
p2 <- make_plot_sim(
  for_lambda = 2**10,
  npoint = 500,
  nsamp = 3
)+xlim(0,10) + ylim(0,10) 


ptot <- egg::ggarrange(p1,p2, nrow=1)
ggsave(plot = ptot, "sim fitted poiss_n3.pdf", width = 4, height = 2)






################################################################################
######################### what is the variance if double sampling? #############
################################################################################


###### Variance decreases in bottleneck   -- normal distribution example
allsamp <- c()
for(i in 1:100){
  dist1 <- rnorm(100000)
  dist2 <- sample(dist1, size=5, replace = TRUE)
  dist3 <- sample(dist2, size=10000, replace = TRUE)
  allsamp <- c(allsamp, sd(dist3))
}
mean(allsamp)


###### Variance decreases in bottleneck   -- poisson distribution example
allsamp <- c()
for(i in 1:500){
  dist1 <- rpois(100000, lambda = 10)
  dist2 <- sample(dist1, size=50, replace = TRUE)
  dist3 <- sample(dist2, size=100000, replace = TRUE)
  allsamp <- c(allsamp, var(dist3))
}
mean(allsamp)

var(dist1)
var(dist2)
var(dist3)



################################################################################
######### Fig S4 - ROC curves ##################################################
################################################################################


#use essential gene list from https://pubmed.ncbi.nlm.nih.gov/30971826/   
depmap_core <- read.csv("depmap_core.csv",sep="\t")
depmap_core <- unique(depmap_core)
depmap_core <- depmap_core[!duplicated(depmap_core$genesym),]  #only MAK16 is there twice
rownames(depmap_core) <- depmap_core$genesym


### Supplemental figure, essentially spread
ggplot(depmap_core, aes(essentiality)) + 
  annotate("rect", xmin = 95, xmax = 100, ymin = 0, ymax = 2500,fill = "lightblue") +
  geom_histogram(breaks = 0:100) + 
  theme_bw() 
ggsave("newout/roc_essentially_spread.svg", width = 5, height = 4)


## Make a ROC from mageck test
make_roc_for_screen_magecktest <- function(dat,withname){
  dat <- dat[dat$id %in% depmap_core$genesym,]
  dim(dat) #note, about half the genes will go missing; but representative set
  dat$is_essential <- depmap_core[dat$id,]$essentiality>0.95
  roc_score <- pROC::roc(dat$is_essential+0, 1-dat$neg.fdr)  
  tc <- pROC::coords(roc_score)
  tc$name <- withname
  tc
}

## Make a ROC from mageck MLE
make_roc_for_screen_mageckMLE <- function(dat,withname){
  print(head(dat))
  dat <- dat[dat$Gene %in% depmap_core$genesym,]
  print(dim(dat))
  #dim(dat) #note, about half the genes will go missing; but representative set
  dat$is_essential <- depmap_core[dat$Gene,]$essentiality>0.95
  #print(head(dat))
  roc_score <- pROC::roc(dat$is_essential+0, dat$eff.z)
  tc <- pROC::coords(roc_score)
  tc$name <- withname
  tc
}



## Plots overlaids ROC from mageck test
make_roc_comp_magecktest <- function(numread){
  toplot <- rbind(
    make_roc_for_screen_magecktest(
      read.csv(file.path(paste0("/corgi/otherdataset/crispr_padlock/original/for_mageck/14d_mip_nodedup#",numread,"#0.gene_summary.txt")),sep="\t"),
      "MIP no dedup"),
    make_roc_for_screen_magecktest(
      read.csv(file.path(paste0("/corgi/otherdataset/crispr_padlock/original/for_mageck/14d_mip_dedup#",numread,"#0.gene_summary.txt")),sep="\t"),
      "MIP dedup"),
    make_roc_for_screen_magecktest(
      read.csv(file.path(paste0("/corgi/otherdataset/crispr_padlock/original/for_mageck/14d_pcr#",numread,"#0.gene_summary.txt")),sep="\t"),
      "PCR")
  )
  ggplot(toplot, aes(specificity, sensitivity, color=name)) + 
    annotate("segment", x = 1, xend = 0, y = 0, yend = 1, colour = "black") +
    geom_path() +
    theme_bw() + 
    scale_x_reverse() +
    theme(legend.title=element_blank())
}
p1 <- make_roc_comp_magecktest("10000") + ggtitle("10k")
p2 <- make_roc_comp_magecktest("100000") + ggtitle("100k")
p3 <- make_roc_comp_magecktest("1000000") + ggtitle("1M")
p4 <- make_roc_comp_magecktest("10000000") + ggtitle("10M")
ptot <- egg::ggarrange(p1,p2,p3,p4, nrow=1)
ptot
ggsave(plot = ptot, "newout/roc_test.svg", width = 20, height = 4)



## Plots overlaids ROC from mageck MLE
make_roc_comp_mageckMLE <- function(numread){
  toplot <- rbind(
    make_roc_for_screen_mageckMLE(
      read.csv(file.path(paste0("/corgi/otherdataset/crispr_padlock/original/for_mageck/14d_mip_nodedup#",numread,"#0.mageckmle.gene_summary.txt")),sep="\t"),
      "MIP no dedup"),
    make_roc_for_screen_mageckMLE(
      read.csv(file.path(paste0("/corgi/otherdataset/crispr_padlock/original/for_mageck/14d_mip_dedup#",numread,"#0.mageckmle.gene_summary.txt")),sep="\t"),
      "MIP dedup"),
    make_roc_for_screen_mageckMLE(
      read.csv(file.path(paste0("/corgi/otherdataset/crispr_padlock/original/for_mageck/14d_pcr#",numread,"#0.mageckmle.gene_summary.txt")),sep="\t"),
      "PCR")
  )
  ggplot(toplot, aes(specificity, sensitivity, color=name)) + 
    annotate("segment", x = 1, xend = 0, y = 0, yend = 1, colour = "black") +
    geom_path() +
    theme_bw() + 
    scale_x_reverse() +
    theme(legend.title=element_blank())
}
p1 <- make_roc_comp_mageckMLE("10000") + ggtitle("10k")
p2 <- make_roc_comp_mageckMLE("100000") + ggtitle("100k")
p3 <- make_roc_comp_mageckMLE("1000000") + ggtitle("1M")
p4 <- make_roc_comp_mageckMLE("10000000") + ggtitle("10M")
ptot <- egg::ggarrange(p1,p2,p3,p4, nrow=1)
ptot
ggsave(plot = ptot, "newout/roc_mle.svg", width = 20, height = 4)









################################################################################
######### Picking of candidates by z-score ##################################### 
################################################################################


#In the Brunello kinome screen, CRISPR-MIP identified ~75 hits and PCR identified 2 hits.
#‘Variance rather than effect size’ was proposed to explain the lack of hits identified using the PCR method.
#MAGeCK hit calling is based on fairly conservative distributional assumptions and gene-level p values are determined by least probable perturbation.
#It would be useful to repeat hit picking with an algorithm less influenced by sgRNA rank (e.g., z-scoring).

#Compare 

dat <- read.csv("/corgi/otherdataset/crispr_padlock/original/for_mageck/14d_pcr#allreads#0.gene_summary.txt" ,sep="\t")

ggplot(dat,aes(neg.fdr, neg.score)) + geom_point()
### is this z score?? cannot believe




dat <- read.csv("/corgi/otherdataset/crispr_padlock/original/for_mageck/14d_pcr#allreads#0.sgrna_summary.txt" ,sep="\t")
dat

############ MLE analysis
dat <- read.csv("/corgi/otherdataset/crispr_padlock/original/for_mageck/14d_pcr#allreads#0.mageckmle.gene_summary.txt" ,sep="\t")
dat

ggplot(dat, aes(eff.z, eff.fdr)) + geom_point()  #went into response letter

####### for MLE command ... so it only works for MLE?
#HL60|p-value	The raw p-value (using permutation) of this gene
#HL60|fdr	The false discovery rate of this gene
#HL60|z	The z-score associated with Wald test













