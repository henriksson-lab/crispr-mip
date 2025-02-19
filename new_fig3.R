 
if(FALSE){
  install.packages("xlsx")
  install.packages("enrichR")
}

library("enrichR")



################################################################################
################## Fig 3A: comparison of mageck TEST/MLE #######################   upper panel
################################################################################

nudge_x <- 100

mageck_mip <- read.csv("/corgi/otherdataset/crispr_padlock/for_mageck/14d_mip_dedup#allreads#0.gene_summary.txt", sep="\t")
mageck_pcr <- read.csv("/corgi/otherdataset/crispr_padlock/for_mageck/14d_pcr#allreads#0.gene_summary.txt", sep="\t")

df <- merge(
  data.frame(
    rank_mip=mageck_mip$neg.rank,
    id=mageck_mip$id,
    fdr_mip=mageck_mip$neg.fdr,
    p_mip=mageck_mip$neg.p.value,
    score_mip=mageck_mip$neg.score
  ),
  data.frame(
    rank_pcr=mageck_pcr$neg.rank,
    id=mageck_pcr$id,
    fdr_pcr=mageck_pcr$neg.fdr,
    p_pcr=mageck_pcr$neg.p.value,
    score_pcr=mageck_pcr$neg.score
  )
)


df$cat <- ""
#df$cat[df$p_mip<0.05 | df$p_pcr<0.05] <- "either"
#df$cat[df$p_mip<0.05 & df$p_pcr<0.05] <- "common"
#df$cat[df$fdr_mip<0.05 | df$fdr_pcr<0.05] <- "either"
#df$cat[df$fdr_mip<0.05 & df$fdr_pcr<0.05] <- "common"
#df$cat[df$fdr_mip<0.05 | df$fdr_pcr<0.05] <- "either"
df$cat[df$fdr_mip<0.05] <- "mip"
df$cat[df$fdr_pcr<0.05] <- "pcr"
df$cat[df$fdr_mip<0.05 & df$fdr_pcr<0.05] <- "common"
df$cat <- factor(df$cat, levels = c("","pcr","mip","common"))
table(df$cat)  ################################################################# all in common

int_genes <- c("CDC7","CDK1","PLK1","DOLK","WEE1","GUK1","AURKB")
#int_genes <- c("CDC7","CDK1","PLK1","DOLK","WEE1","GUK1","AURKB",   "STK32A","FGFRL1")


p1 <- ggplot(df, aes(rank_mip, -log10(score_mip),color=cat, label=id)) + 
  geom_point() +
  theme_bw() +
  scale_colour_manual(values=c("gray","darkgreen","darkblue","red"), breaks=levels(df$cat))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())    +
  ggrepel::geom_text_repel(data=df[df$id %in% int_genes,],nudge_x = nudge_x, max.overlaps = 100)+
  xlab("Gene rank")+
  ylab("-log10(RRA score)")+
  ylim(0,7)+
  theme(legend.position="none")
p1

p2 <- ggplot(df, aes(rank_pcr, -log10(score_pcr),color=cat, label=id)) + 
  geom_point() +
  theme_bw() +
  scale_colour_manual(values=c("gray","darkgreen","darkblue","red"), breaks=levels(df$cat))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())    +
  ggrepel::geom_text_repel(data=df[df$id %in% int_genes,],nudge_x = nudge_x, max.overlaps = 100)+
  xlab("Gene rank")+
  ylab("-log10(RRA score)")+ 
  ylim(0,7)+
  theme(legend.position="none")
p2

#ptot <- egg::ggarrange(p1,p2, nrow=1)



############ number of hits
sum(df$fdr_mip<0.05) #68  MIP test
sum(df$fdr_pcr<0.05) #2   PCR test

################################################################################
################## Fig 3b: GO of mageck TEST/MLE ###############################
################################################################################

################## TEST

mageck_mip <- read.csv("/corgi/otherdataset/crispr_padlock/for_mageck/14d_mip_dedup#allreads#0.gene_summary.txt", sep="\t")
mageck_pcr <- read.csv("/corgi/otherdataset/crispr_padlock/for_mageck/14d_pcr#allreads#0.gene_summary.txt", sep="\t")

setEnrichrSite("Enrichr") # Human genes
dbs <- c("GO_Biological_Process_2023")

enriched_mip <- enrichr(
  genes = mageck_mip$id[mageck_mip$neg.p.value<0.05], 
  background =  mageck_mip$id, 
  dbs
)$GO_Biological_Process_2023

enriched_pcr <- enrichr(
  genes = mageck_pcr$id[mageck_pcr$neg.p.value<0.05], 
  background =  mageck_pcr$id, 
  dbs
)$GO_Biological_Process_2023

go_test_enriched_mip <- enriched_mip
go_test_enriched_pcr <- enriched_pcr

u_enriched_pcr <- enriched_pcr#[!duplicated(enriched_pcr$Genes),]
u_enriched_mip <- enriched_mip#[!duplicated(enriched_mip$Genes),]

toplot <- merge(all=TRUE,
                data.frame(
                  Term=u_enriched_pcr$Term, 
                  p_pcr=u_enriched_pcr$Adjusted.P.value
                ),
                data.frame(
                  Term=u_enriched_mip$Term, 
                  p_mip=u_enriched_mip$Adjusted.P.value
                )
)
toplot$p_pcr[is.na(toplot$p_pcr)] <- 1
toplot$p_mip[is.na(toplot$p_mip)] <- 1

toplot$p_min <- pmin(toplot$p_pcr,toplot$p_mip)
toplot <- toplot[order(toplot$p_min),]
ggplot(toplot[1:30,], aes(Term, -log10(p_mip))) +
  geom_bar(stat="identity") + 
  coord_flip()

ggplot(toplot[1:30,], aes(Term, -log10(p_pcr))) +
  geom_bar(stat="identity") + 
  coord_flip()


toplot_test <- toplot
pscatter_test <- ggplot(toplot_test, aes(p_pcr, p_mip)) + 
  geom_point() + 
  theme_bw() +
  xlab("p PCR") + ylab("p MIP") + xlim(0,1) + ylim(0,1)




################## MLE
mageck_mip <- read.csv("/corgi/otherdataset/crispr_padlock/for_mageck/14d_mip_dedup#allreads#0.mageckmle.gene_summary.txt", sep="\t")
mageck_pcr <- read.csv("/corgi/otherdataset/crispr_padlock/for_mageck/14d_pcr#allreads#0.mageckmle.gene_summary.txt", sep="\t")

setEnrichrSite("Enrichr") # Human genes
#listEnrichrDbs()
dbs <- c("GO_Biological_Process_2023")

enriched_mip <- enrichr(
  genes = mageck_mip$Gene[mageck_mip$eff.p.value<0.05], 
  background =  mageck_mip$Gene, 
  dbs
)$GO_Biological_Process_2023

enriched_pcr <- enrichr(
  genes = mageck_pcr$Gene[mageck_pcr$eff.p.value<0.05], 
  background =  mageck_pcr$Gene, 
  dbs
)$GO_Biological_Process_2023

go_mle_enriched_mip <- enriched_mip
go_mle_enriched_pcr <- enriched_pcr

u_enriched_pcr <- enriched_pcr#[!duplicated(enriched_pcr$Genes),]
u_enriched_mip <- enriched_mip#[!duplicated(enriched_mip$Genes),]

toplot <- merge(all=TRUE,
                data.frame(
                  Term=u_enriched_pcr$Term, 
                  p_pcr=u_enriched_pcr$Adjusted.P.value
                ),
                data.frame(
                  Term=u_enriched_mip$Term, 
                  p_mip=u_enriched_mip$Adjusted.P.value
                )
)
toplot$p_pcr[is.na(toplot$p_pcr)] <- 1
toplot$p_mip[is.na(toplot$p_mip)] <- 1

toplot$p_min <- pmin(toplot$p_pcr,toplot$p_mip)
toplot <- toplot[order(toplot$p_min),]
ggplot(toplot[1:30,], aes(Term, -log10(p_mip))) +
  geom_bar(stat="identity") + 
  coord_flip() 

ggplot(toplot[1:30,], aes(Term, -log10(p_pcr))) +
  geom_bar(stat="identity") + 
  coord_flip() 
#egg::ggarrange(p1,p2,p3,p4,nrow=2)


toplot_mle <- toplot
pscatter_mle <- ggplot(toplot, aes(p_pcr, p_mip)) + 
  geom_point() + 
  theme_bw() +
  xlab("p PCR") + ylab("p MIP") + xlim(0,1) + ylim(0,1)
  #scale_x_log10() + scale_y_log10()

ptot <- egg::ggarrange(pscatter_test,pscatter_mle,nrow=1)
ptot
ggsave(plot = ptot, "newout/scatter_go_fig3.svg", width = 6, height = 3)


################################################################################
################## Fig 3A: comparison of mageck TEST/MLE #######################   lower panel
################################################################################


mageck_mip <- read.csv("/corgi/otherdataset/crispr_padlock/for_mageck/14d_mip_dedup#allreads#0.mageckmle.gene_summary.txt", sep="\t")
mageck_pcr <- read.csv("/corgi/otherdataset/crispr_padlock/for_mageck/14d_pcr#allreads#0.mageckmle.gene_summary.txt", sep="\t")

df <- merge(
  data.frame(
    fdr_mip=mageck_mip$eff.fdr,
    rank_mip=rank(mageck_mip$eff.beta),
    beta_mip=mageck_mip$eff.beta,
    id=mageck_mip$Gene
  ),
  data.frame(
    fdr_pcr=mageck_pcr$eff.fdr,
    rank_pcr=rank(mageck_pcr$eff.beta), 
    beta_pcr=mageck_pcr$eff.beta,
    id=mageck_pcr$Gene
  )
)

df$cat <- ""
df$cat[df$fdr_mip<0.05] <- "mip"
df$cat[df$fdr_pcr<0.05] <- "pcr"
df$cat[df$fdr_mip<0.05 & df$fdr_pcr<0.05] <- "common"
df$cat <- factor(df$cat, levels = c("","pcr","mip","common"))

int_genes <- c("CDC7","CDK1","PLK1","DOLK","WEE1","GUK1","AURKB")

df[df$id %in% int_genes,]



######## MIP plot
p3 <- ggplot(df, aes(rank_mip, beta_mip,color=cat, label=id)) + 
  geom_point() +
  theme_bw() +
  scale_colour_manual(values=c("gray","darkgreen","darkblue","red"), breaks=levels(df$cat))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())    +
  ggrepel::geom_text_repel(data=df[df$id %in% int_genes,],nudge_x = nudge_x, max.overlaps = 100)+
  xlab("Gene rank")+
  ylab("Beta score")+
  ylim(-2.1,0.7)+
  theme(legend.position="none")
p3


######## PCR plot
p4 <- ggplot(df, aes(rank_pcr, beta_pcr,color=cat, label=id)) + 
  geom_point() +
  theme_bw() +
  scale_colour_manual(values=c("gray","darkgreen","darkblue","red"), breaks=levels(df$cat))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())    +
  ggrepel::geom_text_repel(data=df[df$id %in% int_genes,],nudge_x = nudge_x, max.overlaps = 100)+
  xlab("Gene rank")+
  ylab("Beta score")+
  ylim(-2.1,0.7)+
  theme(legend.position="none")
p4

ptot <- egg::ggarrange(p1,p2,p3,p4, nrow=2)
ptot
ggsave(plot=ptot, "out/fig3a.svg", width = 6, height = 6)


sum(df$fdr_mip<0.05) #63  MIP test (mle)
sum(df$fdr_pcr<0.05) #20  PCR test (mle)




################################################################################
################## Fig 3C: comparison of mageck TEST/MLE #######################
################################################################################

int_genes <- c("AURKB", "DOLK", "GUK1", "WEE1")

######### test
######### test
######### test

mageck_mip <- read.csv("/corgi/otherdataset/crispr_padlock/for_mageck/14d_mip_dedup#allreads#0.gene_summary.txt", sep="\t")
mageck_pcr <- read.csv("/corgi/otherdataset/crispr_padlock/for_mageck/14d_pcr#allreads#0.gene_summary.txt", sep="\t")

df <- merge(
  data.frame(
    rank_mip=mageck_mip$neg.rank,
    id=mageck_mip$id,
    fdr_mip=mageck_mip$neg.fdr
  ),
  data.frame(
    rank_pcr=mageck_pcr$neg.rank,
    id=mageck_pcr$id,
    fdr_pcr=mageck_pcr$neg.fdr
  )
)
df$score <- df$fdr_mip*df$fdr_pcr

df$cat <- ""
df$cat[df$id %in% int_genes] <- "fo"

df <- df[order(df$cat),]
p1 <- ggplot(df, aes(rank_mip, rank_pcr, color=cat)) + 
  geom_point() +
  scale_colour_manual(values=c("lightblue","red"))+
  xlab("Rank RRA (CRISPR-MIP)") +
  ylab("Rank RRA (CRISPR-PCR)") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())  +
  theme(legend.position = "none")
p1

cor.test(df$rank_mip, df$rank_pcr) #0.79


######### MLE
######### MLE
######### MLE

mageck_mip <- read.csv("/corgi/otherdataset/crispr_padlock/for_mageck/14d_mip_dedup#allreads#0.mageckmle.gene_summary.txt", sep="\t")
mageck_pcr <- read.csv("/corgi/otherdataset/crispr_padlock/for_mageck/14d_pcr#allreads#0.mageckmle.gene_summary.txt", sep="\t")

df <- merge(
  data.frame(
    #z_mip=mageck_mip$eff.z,
    rank_mip=rank(mageck_mip$eff.z),
    beta_mip=mageck_mip$eff.beta,
    id=mageck_mip$Gene
  ),
  data.frame(
    #z_pcr=mageck_pcr$eff.z,
    rank_pcr=rank(mageck_pcr$eff.z), 
    beta_pcr=mageck_pcr$eff.beta,
    id=mageck_pcr$Gene
  )
)

df$cat <- ""
df$cat[df$id %in% int_genes] <- "fo"
df <- df[order(df$cat),]

p2 <- ggplot(df, aes(rank_mip, rank_pcr, color=cat)) + 
  geom_point() +
  scale_colour_manual(values=c("lightblue","red"))+
  xlab("Rank z-score (CRISPR-MIP)") +
  ylab("Rank z-score (CRISPR-PCR)") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")
p2

cor.test(df$rank_mip, df$rank_pcr) #0.83


ptot <- egg::ggarrange(p1,p2)
ptot
ggsave(plot = ptot, "out/fig3c.pdf", width = 4, height = 6) 




################################################################################
################## Fig 3D: sgRNA abundances ####################################
################################################################################


list_gene <- c("AURKB","CDK1","DOLK","GUK1","WEE1")

sgrna_abundance  <- read.csv("/corgi/otherdataset/crispr_padlock/for_mageck/14d_mip_dedup#allreads#0.sgrna_summary.txt",sep="\t")

df <- rbind(
  data.frame(
    sgrna=sgrna_abundance$sgrna,
    gene=sgrna_abundance$Gene,
    cnt=sgrna_abundance$control_mean,
    time="d0" 
  ),
  data.frame(
    sgrna=sgrna_abundance$sgrna,
    gene=sgrna_abundance$Gene,
    cnt=sgrna_abundance$treat_mean,
    time="d14" 
  )
)
df$time <- factor(df$time, levels = c("d14","d0"))

all_plots <- list()

max_scale <- max(df[df$gene %in% list_gene,]$cnt)
for(curg in list_gene){
  all_plots[[paste("mip", curg)]] <- ggplot(
    df[df$gene==curg,], 
    aes(sgrna, time, fill=cnt)) + 
    geom_tile(color="black") + 
    scale_fill_gradientn(colours = c("blue","white","red"), limits=c(0,max_scale)) +
    theme_bw()+ theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_blank())+  #element_line(colour = "black")
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    ggtitle(curg)
}






sgrna_abundance  <- read.csv("/corgi/otherdataset/crispr_padlock/for_mageck/14d_pcr#allreads#0.sgrna_summary.txt",sep="\t")

df <- rbind(
  data.frame(
    sgrna=sgrna_abundance$sgrna,
    gene=sgrna_abundance$Gene,
    cnt=sgrna_abundance$control_mean,
    time="d0" 
  ),
  data.frame(
    sgrna=sgrna_abundance$sgrna,
    gene=sgrna_abundance$Gene,
    cnt=sgrna_abundance$treat_mean,
    time="d14" 
  )
)
df$time <- factor(df$time, levels = c("d14","d0"))

max_scale <- max(df[df$gene %in% list_gene,]$cnt)
for(curg in list_gene){
  all_plots[[paste("pcr", curg)]] <- ggplot(
    df[df$gene==curg,], 
    aes(sgrna, time, fill=cnt)) + 
    geom_tile(color="black") + 
    scale_fill_gradientn(colours = c("blue","white","red"), limits=c(0,max_scale)) +
    theme_bw()+ theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_blank())+  
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    ggtitle(curg)
}


ptot <- egg::ggarrange(plots=all_plots, nrow=2)
ptot
ggsave(filename = "newout/fig3d.svg", plot=ptot, width = 18, height = 3.5)
  

################################################################################
################## Create a joint XLS sheet with data ##########################
################################################################################

library("xlsx")

dir_mageck <- "/corgi/otherdataset/crispr_padlock/for_mageck"
mageck_mip <- read.csv("/corgi/otherdataset/crispr_padlock/for_mageck/14d_mip_dedup#allreads#0.mageckmle.gene_summary.txt", sep="\t")
mageck_pcr <- read.csv("/corgi/otherdataset/crispr_padlock/for_mageck/14d_pcr#allreads#0.mageckmle.gene_summary.txt", sep="\t")


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

out_xls <- "/corgi/otherdataset/crispr_padlock/out_sup_table1.xlsx"

for(i in 1:length(list_cond)){
  dat <- read.csv(file.path(dir_mageck, list_fname[i]),sep="\t", as.is = TRUE, check.names = FALSE)
  if("neg|fdr" %in% colnames(dat)){
    numsig <- sum(dat$`neg|fdr`<0.05)
  } else {
    numsig <- sum(dat$`eff|fdr`<0.05)
    dat <- dat[order(dat$`eff|fdr`<0.05),]
  }
  print(paste(list_cond[i], "---", list_fname[i], "--- num sig:", numsig))
  write.xlsx(dat, file = out_xls,sheetName = list_cond[i], append = i!=1, row.names = FALSE)  
}


write.xlsx(go_mle_enriched_mip, file = out_xls, sheetName = "GO MLE MIP", append = TRUE, row.names = FALSE)  
write.xlsx(go_mle_enriched_pcr, file = out_xls, sheetName = "GO MLE PCR", append = TRUE, row.names = FALSE)  
write.xlsx(go_test_enriched_mip, file = out_xls, sheetName = "GO test MIP", append = TRUE, row.names = FALSE)
write.xlsx(go_test_enriched_pcr, file = out_xls, sheetName = "GO test PCR", append = TRUE, row.names = FALSE)



