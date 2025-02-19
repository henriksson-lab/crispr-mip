library(stringr)
#library(DESeq2)
library(ggplot2)
library(dplyr)



################################################################################
########### Helper function ####################################################
################################################################################

gcproc <- function(grna){
  (stringr::str_count(grna, "G") + stringr::str_count(grna, "C")) / str_length(grna)
}


################################################################################
############ Fig4a: sgRNA abundances, PCR vs padlock, for one hit ##############
################################################################################


plotSgrnaAbundance <- function(dat){
  df <- rbind(
    data.frame(
      sgrna=dat$sgrna,
      gene=dat$Gene,
      mean=dat$control_mean,
      sample="lib"
    ),
    data.frame(
      sgrna=dat$sgrna,
      gene=dat$Gene,
      mean=dat$treat_mean,
      sample="treatment"
    )
  )
  ggplot(df, aes(sample,mean, group=sgrna)) + geom_line() + theme_bw()
}



dat <- read.csv("/corgi/otherdataset/crispr_padlock/for_mageck/14d_mip_dedup#4000000#0.sgrna_summary.txt",sep="\t")
dat$control_mean <- dat$control_mean/sum(dat$control_mean)*1e6
dat$treat_mean <- dat$treat_mean/sum(dat$treat_mean)*1e6
p1 <- plotSgrnaAbundance(dat[dat$Gene=="CDK1",]) + xlab("MIP dedup") + ylab("Rel. abundance") + ylim(0,400)

dat <- read.csv("/corgi/otherdataset/crispr_padlock/for_mageck/14d_pcr#4000000#0.sgrna_summary.txt",sep="\t")
dat$control_mean <- dat$control_mean/sum(dat$control_mean)*1e6
dat$treat_mean <- dat$treat_mean/sum(dat$treat_mean)*1e6
p2 <- plotSgrnaAbundance(dat[dat$Gene=="CDK1",]) + xlab("PCR") + ylab("Rel. abundance") + ylim(0,400)

ptot <- egg::ggarrange(p1,p2, nrow = 1)
ptot 

ggsave("out/fig4a.svg", ptot, width = 5, height = 2.5)




################################################################################
#################### Get metadata for full-depth samples #######################
################################################################################



#update later
geneinfo <- read.csv("/corgi/otherdataset/crispr_padlock/merged/lib.csv", sep="\t")
geneinfo <- data.frame(
  gene=geneinfo$Target.Gene.Symbol,
  grna=geneinfo$sgRNA.Target.Sequence
)

samplemeta_all <- read.csv("/husky/martin/CRISPR_padlock_all_fastq/samplemeta.csv",sep="\t")
samplemeta_pcr   <- samplemeta_all[samplemeta_all$protocol=="PCR",]
samplemeta_padlock_dedup <- samplemeta_all[samplemeta_all$protocol=="padlock",]
samplemeta_padlock_nodedup <- samplemeta_all[samplemeta_all$protocol=="padlock",]

samplemeta_pcr$filename <- paste0(samplemeta_pcr$filename, ".csv")
samplemeta_padlock_dedup$filename <- paste0(samplemeta_padlock_dedup$filename, ".dedup.csv")
samplemeta_padlock_nodedup$filename <- paste0(samplemeta_padlock_nodedup$filename, ".nodedup.csv")

samplemeta_pcr$dedup <- FALSE
samplemeta_padlock_dedup$dedup <- TRUE
samplemeta_padlock_nodedup$dedup <- FALSE

samplemeta <- rbind(
  samplemeta_pcr,
  samplemeta_padlock_dedup,
  samplemeta_padlock_nodedup
)

################################################################################
############# Fig 4b new : why are padlock and pcr libs different? #############
################################################################################


############# Get counts for full-depth sample  (ignore non-dedup)
alldat <- list()
basedir <- "/corgi/otherdataset/crispr_padlock/subsamp/"
for(f in list.files(basedir)){
  if(str_ends(f, fixed("_dedup#allreads#0")) || str_ends(f, fixed("_pcr#allreads#0"))) {   
    dat <- read.csv(file.path(basedir,f))
    
    f <- str_remove(f, fixed("_dedup#allreads#0"))
    f <- str_remove(f, fixed("_pcr#allreads#0"))
    print(f)

    dat$filename <- f
    alldat[[f]] <- dat
  }
}
dat <- do.call(rbind, alldat)

samplemeta <- read.csv("/corgi/otherdataset/crispr_padlock/martin_samplemeta.csv",sep="\t")
dat <- merge(dat, samplemeta)

df <- merge(
  sqldf("select avg(cnt) as cnt_mip, grna from dat where time='14d' and protocol='padlock' group by grna"),
  sqldf("select avg(cnt) as cnt_pcr, grna from dat where time='14d' and protocol='PCR' group by grna"),
)
df$gc <- gcproc(df$grna)  

p1 <- ggplot(df, aes(cnt_pcr, cnt_mip, color=gc)) + 
  geom_point() + 
  scale_x_log10()+
  scale_y_log10()+
  theme_bw()+
  xlab("PCR abundance") + 
  ylab("MIP dedup abundance")+
  scale_color_gradientn(colors = c("red","gray","blue"),values=c(0,0.5,1))+
  geom_abline(intercept = -1.8, slope = 1)


df <- merge(
  sqldf("select avg(cnt) as cnt_mip, grna from dat where time='lib' and protocol='padlock' group by grna"),
  sqldf("select avg(cnt) as cnt_pcr, grna from dat where time='lib' and protocol='PCR' group by grna"),
)
df$gc <- gcproc(df$grna)  

p2 <- ggplot(df, aes(cnt_pcr, cnt_mip, color=gc)) + 
  geom_point() + 
  scale_x_log10()+
  scale_y_log10()+
  theme_bw()+
  xlab("PCR abundance") + 
  ylab("MIP dedup abundance")+
  scale_color_gradientn(colors = c("red","gray","blue"),values=c(0,0.5,1))+
  geom_abline(intercept = -0.05, slope = 1)


ptot <- egg::ggarrange(p1,p2, nrow=1)
ggsave(plot = ptot, "out/fig4b_new.svg", width = 7, height = 3)



################################################################################
################## Fig 4b why are padlock and pcr libs different? ##############
################################################################################


df <- data.frame(
  pcr=rowMeans(cnt[,samplemeta$time=="lib" & samplemeta$protocol=="PCR"]),
  padlock=rowMeans(cnt[,samplemeta$time=="lib" & samplemeta$protocol=="padlock"]),
  grna=rownames(cnt)
)
df$gc <- gcproc(df$grna)  
ggplot(df, aes(pcr, padlock, label=grna, color=gc)) + 
  geom_point() + 
  scale_x_log10() + scale_y_log10() +
  xlab("PCR abundance") + ylab("MIP dedup abundance")+
  scale_color_gradientn(colors = c("red","gray","blue"),values=c(0,0.5,1))+
  geom_abline(intercept = 0.6, slope = 1)
ggsave("out/fig4b.svg", width = 5, height = 3)



################################################################################
################## fig1e efficiency vs padlock amount ##########################
################################################################################


#Only get subset of samples that are relevant
samplemeta_conc <- samplemeta_padlock_dedup[
    str_starts(samplemeta_padlock_dedup$filename,"padlock_A") | str_starts(samplemeta_padlock_dedup$filename,"padlock_B") | str_starts(samplemeta_padlock_dedup$filename,"padlock_C") |
       samplemeta_padlock_dedup$library=="GFPg1"
    ,]

#Confirmed number of cells to be correct. 150k for gfp, 1M for other efficiency test
samplemeta_conc$conc <- as.double(str_split_i(str_split_i(samplemeta_conc$samplename,"gDNA_",2),"_",1))


alldf <- NULL
for(i in 1:nrow(samplemeta_conc)){
  
  ftemp <- samplemeta_conc$filename[i]
  ftemp <- str_remove(ftemp, fixed(".dedup.csv"))
  ftemp <- paste0(ftemp, "_dedup#allreads#0")
  cur_file <- file.path("/corgi/otherdataset/crispr_padlock/subsamp",ftemp) 
  dat_dedup <- read.csv(cur_file)
  nrow(dat_dedup)
  
  ftemp <- samplemeta_conc$filename[i]
  ftemp <- str_remove(ftemp, fixed(".dedup.csv"))
  ftemp <- paste0(ftemp, "_nodedup#allreads#0")
  cur_file <- file.path("/corgi/otherdataset/crispr_padlock/subsamp",ftemp)
  dat_nodedup <- read.csv(cur_file)
  

  num_reads <- sum(dat_nodedup$cnt)
  cnt_dedup <- sum(dat_dedup$cnt)
  num_cells <- samplemeta_conc$cells[i]
  
  should_obs <- length(unique(sample(x = 1:num_cells, size = num_reads, replace = TRUE)))
  
  
  alldf <- rbind(
    alldf,
    data.frame(
      filename=samplemeta_conc$filename[i], 
      cnt_dedup=cnt_dedup,                     
      num_reads=num_reads,
      should_obs=should_obs
    )
  )
  
}
alldf$ratio <- alldf$cnt_dedup/alldf$should_obs
alldf <- merge(samplemeta_conc,alldf)

alldf$is_digested <- !str_detect(alldf$treatment, "non-digested")
alldf$conc[alldf$library=="GFPg1"] <- 0.02

alldf

ggplot(alldf, aes(conc, ratio*100,color=is_digested, shape=library)) +   
  geom_point() + 
  scale_x_log10() + 
  theme_bw() +
  ylab("% captured") + 
  xlab("pM padlock")
ggsave("newout/fig1e.svg", width = 4, height = 2.5)



################################################################################
########### Fig 4c: Compare logFC between MIP and PCR, some genes ##############
################################################################################

### Compute GC average per gene
gene2grna <- read.csv("/corgi/otherdataset/crispr_padlock/for_mageck/14d_mip_dedup#4000000#0.sgrna_summary.txt",sep="\t")
gene2grna$gc <- gcproc(gene2grna$sgrna)
gene_gc <- sqldf::sqldf("select sgrna, Gene as gene, avg(gc) as gc from gene2grna group by gene")



dat1 <- read.csv("/corgi/otherdataset/crispr_padlock/for_mageck/14d_mip_dedup#4000000#0.gene_summary.txt",sep="\t")
dat2 <- read.csv("/corgi/otherdataset/crispr_padlock/for_mageck/14d_pcr#4000000#0.gene_summary.txt",sep="\t")
df <- merge(
  data.frame(
    gene=dat1$id,
    log_fc_mipdedup=dat1$pos.lfc
  ),
  data.frame(
    gene=dat2$id,
    log_fc_pcr=dat2$pos.lfc
  )
)
df <- merge(df, gene_gc) ####### missing!
df$is_aurkb <- df$gene=="AURKB"
df <- df[order(df$is_aurkb),]
ggplot(df, aes(log_fc_mipdedup, log_fc_pcr, color=gc)) + 
  geom_point() + title("14d 4M reads") + 
  xlab("log(FC) MIPdedup, per gene") + ylab("log(FC) PCR, per gene") +
  scale_color_gradientn(colors = c("red","gray","blue"),values=c(0,0.5,1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.background = element_blank())
ggsave("out/fig4c.svg", width = 3.5, height = 3)




################################################################################
################## Figure 4de - Analysis of Sintov2022 #########################
################################################################################



library(DESeq2)
library(stringr)
library(ggplot2)

dat <- read.csv("sintov2022/GSM6008413_counts-JD_GPP2373_2901507_Sintov_20210224.txt.gz",sep="\t")
cnt <- dat[,-c(1:3)]
cnt <- cnt[,colnames(cnt)!="NTC"]

cond <- data.frame(
  row.names = colnames(cnt)
)

dds <- DESeq2::DESeqDataSetFromMatrix(countData = cnt,
                                      colData = cond,
                                      design= ~1)
dds <- DESeq(dds)
ncnt <- counts(dds, normalized=TRUE)

df <- data.frame(
  x=ncnt[,"Ms..64.Ctrl"],
  y=ncnt[,"CP0043.Library.DNA"]
)
df$fc <- log10(df$x/df$y)
df$gc <- gcproc(dat$Construct.Barcode)
df$gc.jitter <- df$gc + runif(nrow(df), min = -0.5/20, max=0.5/20)

if(FALSE){
  #ggplot(df, aes(fc, gc)) + geom_point() + stat_smooth(method=lm)
  ggplot(df, aes(fc, gc)) + geom_point(mapping = aes(fc, gc.jitter)) + stat_smooth(method=lm)
  ggplot(df, aes(fc, gc)) + geom_bin2d(mapping = aes(fc, gc.jitter),bins=50) + stat_smooth(method=lm)
  
  
  cor.test(df$LFC,df$gc, method = "spearman")
  cor.test(df$LFC[dat$NTC==1],df$gc[dat$NTC==1], method = "spearman")
  
  p1 <- ggplot(df, aes(fc, gc)) + 
    geom_bin2d(bins=50) + stat_smooth(method=lm) + xlab("log10 FC, ctrl/lib, all sgRNA")
  p2 <- ggplot(df[dat$NTC==1,], aes(fc, gc)) + 
    geom_bin2d(bins=50) + stat_smooth(method=lm) + xlab("log10 FC, ctrl/lib, NTC sgRNA")
  ptot <- egg::ggarrange(p1, p2, ncol=2)
  ptot

}


p1 <- ggplot(df,              aes(fc, gc)) + geom_bin2d(mapping = aes(fc, gc.jitter),bins=50) + stat_smooth(method=lm) + xlab("log10 FC, ctrl/lib, all sgRNA") + ylab("GC")
p2 <- ggplot(df[dat$NTC==1,], aes(fc, gc)) + geom_bin2d(mapping = aes(fc, gc.jitter),bins=50) + stat_smooth(method=lm) + xlab("log10 FC, ctrl/lib, NTC sgRNA") + ylab("GC")
ptot <- egg::ggarrange(p1, p2, ncol=2)
ptot
ggsave(plot = ptot, "out/fig4de.svg", width = 8, height = 3)





################################################################################
################################################################################
################################################################################
################## Figure 4 - Depmap analysis ##################################
################################################################################
################################################################################
################################################################################

library(sn)


depmap_grna <- read.csv("/husky/otherdataset/depmap/grna.csv",sep="\t")
head(depmap_grna)
colnames(depmap_grna) <- c("sgrna","Gene","grna_seq","hgnc")
depmap_grna$gc <- gcproc(depmap_grna$grna_seq) + runif(nrow(depmap_grna),-0.5/str_length("CCGAACAAAGTAGCAGTGC"), +0.5/str_length("CCGAACAAAGTAGCAGTGC"))
depmap_grna$gc.jitter <- depmap_grna$gc + runif(nrow(depmap_grna),-0.5/str_length("CCGAACAAAGTAGCAGTGC"), +0.5/str_length("CCGAACAAAGTAGCAGTGC"))


#######################
####################### Average all screens -- Fig 4f
#######################

allf <- list()
for(f in list.files("/husky/otherdataset/depmap/00_raw_counts/V1.0_FirstBatch/magecktest/", pattern = "*.sgrna_summary.txt")){
  print(f)
  dat <- read.csv(file.path("/husky/otherdataset/depmap/00_raw_counts/V1.0_FirstBatch/magecktest/",f),sep="\t")
  dat <- dat[,c("sgrna","Gene","LFC","control_count","treatment_count")]
  dat$f <- f
  allf[[f]] <- dat
}
allfm <- do.call(rbind,allf)

allf_avg <- sqldf::sqldf("select sgrna, Gene, avg(LFC) as LFC, avg(control_count) as control_count, avg(treatment_count) as treatment_count from allfm group by sgrna, Gene")  #possibly better to cast to matrix first instead
allf_avg <- merge(allf_avg, depmap_grna)

my_breaks <- c(1,5,10, 50, 100,500, 1000)

if(FALSE){
  ggplot(allf_avg, aes(LFC,gc-mean(gc))) + 
    geom_bin2d(aes(LFC,gc-mean(gc)), bins = 50) + 
    xlab("Average LFC, DepMap, 356 screens v1.0 set") + 
    geom_smooth(method = "lm", formula = y ~ x+0) + 
    scale_fill_gradient(name = "count", trans = "log",breaks = my_breaks, labels = my_breaks)
  
  cor.test(allf_avg$LFC,allf_avg$gc, method = "spearman")
  
}

#######################
####################### Compare distributions
#######################

library(sn)
mod <- selm(LFC ~ 1, data=allf_avg)
extractSECdistr(mod)

#library(fitdistrplus)
#ca <- fitdistr(allf_avg$LFC, "cauchy")

#######################
####################### Skew normal fit
#######################

if(FALSE){
  thex <- seq(from=-5, to=2, by=0.01)
  plot(density(allf_avg$LFC))
  lines(thex, sn::dsn(seq(from=-5, to=2, by=0.01), xi=0.4711106, omega=1.0980625, alpha= -8.2328289), col="red")
  lines(thex, dcauchy(thex, location=-0.021902958,    scale=0.237180401 ), col="blue")
  
  
}

######################## Fit model
df <- data.frame(
  LFC=allf_avg$LFC,
  gc=allf_avg$gc 
)
df$gc2 <- df$gc**2
df$gc3 <- df$gc**3
df$gc4 <- df$gc**4

mod <- selm(LFC ~ gc + gc2 + gc3 + gc4, data=df)  #, family = "ST")
df$LFC.pred <- predict(mod)
sn_param <- extractSECdistr(mod)

df$lfc <- rsn(nrow(df), xi=sn_param@dp["xi"], omega=sn_param@dp["omega"], alpha=sn_param@dp["alpha"])


######################## Plot fitted model
ggplot(allf_avg, aes(LFC,gc.jitter)) + 
  geom_bin2d(bins = 50) + 
  xlab("Average LFC, DepMap, 356 screens v1.0 set") + 
  geom_smooth(method = "lm", formula = y ~ x) + 
  scale_fill_gradient(name = "count", trans = "log",breaks = my_breaks, labels = my_breaks) +
  geom_point(data = df, mapping=aes(LFC.pred-median(df$LFC.pred), gc))
ggsave("out/fig4f.svg", width=5, height=4)

ggplot(allf_avg, aes(LFC,gc.jitter)) + 
  geom_bin2d(bins = 50) + 
  xlab("Average LFC, DepMap, 356 screens v1.0 set") + 
  scale_fill_gradient(name = "count", trans = "log",breaks = my_breaks, labels = my_breaks) +
  geom_point(data = df, mapping=aes(LFC.pred, gc))
ggsave("out/fig4f_overlay.svg", width=5, height=4)




mod_coef <- as.list(coef(mod))

if(FALSE){
  ggplot(allf_avg, aes(LFC,gc.jitter)) + 
    geom_bin2d(bins = 50) + 
    xlab("Average LFC, DepMap, 356 screens v1.0 set") + 
    geom_smooth(method = "lm", formula = y ~ x) + 
    scale_fill_gradient(name = "count", trans = "log",breaks = my_breaks, labels = my_breaks) +
    geom_point(data = df, mapping=aes(
      mod_coef[[1]] + mod_coef$gc*gc + mod_coef$gc2*gc2 + mod_coef$gc3*gc3 + mod_coef$gc4*gc4 - median(df$LFC.pred), 
      gc))
  
}



#######################
####################### Do correction given GC model
#######################

predLfcCompensation <- function(gc){
  gc2 <- gc**2
  gc3 <- gc**3
  gc4 <- gc**4
  mod_coef[[1]] + mod_coef$gc*gc + mod_coef$gc2*gc2 + mod_coef$gc3*gc3 + mod_coef$gc4*gc4 - median(df$LFC.pred)
}

#LFC is log2

allf_avg$LFC_corrected <- allf_avg$LFC - predLfcCompensation(allf_avg$gc)

ggplot(allf_avg, aes(LFC_corrected,gc.jitter)) + 
  geom_bin2d(bins = 50) + 
  xlab("Average LFC, DepMap, 356 screens v1.0 set -- GC corrected") + 
  #geom_smooth(method = "lm", formula = y ~ x) + 
  scale_fill_gradient(name = "count", trans = "log",breaks = my_breaks, labels = my_breaks) 
ggsave("out/fig4g.svg", width=5, height=4)





#######################
####################### Do correction given GC model --  on counts
#######################

toavg <- read.csv("/husky/otherdataset/depmap/00_raw_counts/V1.0_FirstBatch/magecktest/MIAPac_c904R1.read_count.tsv.gz.sgrna_summary.txt",sep="\t")
toavg <- merge(toavg, depmap_grna)
toavg$LFC_corrected <- toavg$LFC - predLfcCompensation(allf_avg$gc)

gene_avg_lfc <- merge(
  sqldf::sqldf("select Gene, avg(LFC) as LFC_uncorr, count(LFC) as sgrna_cnt from toavg group by Gene"),
  sqldf::sqldf("select Gene, avg(LFC_corrected) as LFC_corr   from toavg group by Gene")
)
gene_avg_lfc$delta <- abs(gene_avg_lfc$LFC_uncorr - gene_avg_lfc$LFC_corr)
gene_avg_lfc$mispredicted <- gene_avg_lfc$delta>0.4
sum(gene_avg_lfc$mispredicted)

ggplot(gene_avg_lfc, aes(LFC_uncorr, LFC_corr, label=Gene)) + geom_point(color="gray") + 
  ggrepel::geom_text_repel(data=gene_avg_lfc[gene_avg_lfc$mispredicted,], color="red", max.overlaps = 1000, segment.colour="black") +
  xlab("LFC uncorrected") + ylab("LFC corrected")
ggsave("out/fig4i.svg", width=5, height=4)


gene_avg_lfc$sgrna_cnt.jitter <- gene_avg_lfc$sgrna_cnt + runif(nrow(gene_avg_lfc), -0.5,0.5)
ggplot(gene_avg_lfc, aes(delta, sgrna_cnt.jitter, label=Gene)) + 
  geom_point(color="black") +
  xlab("Delta log2 FC") + ylab("sgRNA count") 
ggsave("out/fig4j.svg", width=5, height=4)












