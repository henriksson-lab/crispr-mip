
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
################## Figure 1efg + supfig3 - screen hits vs depth ################
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
    numread <- as.integer(str_split_fixed(f,"#",3)[2])
    
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


########## plot with conf interval - d14
grouped <- group_by(out_summary[!is.na(out_summary$numread) & out_summary$numread>min_read_cutoff,], cond, numread)
forgg <- summarise(grouped, mean=mean(numsig), sd=sd(numsig))
p_grna_test_14d <- ggplot(data = forgg[str_detect(forgg$cond,"14"),], aes(x = numread, group = cond)) + 
  geom_line(aes(y = mean, color = cond), size = 1) + 
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = cond), alpha = .2) +
  xlab("#reads") + ylab("MAGeCK TEST #sgRNA (d14)") + scale_x_log10()

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





################### plotting: mageck MLE analysis - genes

out_summary <- NULL
listsum <- list.files(dir_mageck, pattern = "*.mageckmle.gene_summary.txt")
for(f in listsum){
  curcond <- str_split_fixed(f,"#",3)[1]
  numread <- as.integer(str_split_fixed(f,"#",3)[2])
  
  dat <- read.csv(file.path(dir_mageck, f),sep="\t")
  numsig <- sum(dat$eff.p.value<0.05)
  
  df <- data.frame(
    cond=curcond,
    numread=numread,
    numsig=numsig
  )
  out_summary <- rbind(df, out_summary)
}


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
  format_one_panel(p_grna_test_14d)+ylim(0,1000),
  format_one_panel(p_gene_test_14d)+ylim(0,200),
  format_one_panel(p_gene_mle_14d)+ylim(0,200),
  ncol = 3
)
ptot
ggsave(filename = "out/fig1_efg.svg", ptot, width = 6*2, height = 3*3)


