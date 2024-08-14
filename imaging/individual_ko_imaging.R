
fname_input_data <- "individual_ko_imaging.csv"


## Ensure consistent coloring
dat <- read.csv(fname_input_data,sep="\t")
all_gene_names <- unique(dat[1,-(1:4)])
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
gene_colors <- gg_color_hue(length(all_gene_names))
names(gene_colors) <- all_gene_names

gene_colors["NTC"] <- "#AAAAFF"

dat

#################
#################  Plot all samples with normalization
#################


allplot <- list()
all_reps <- c("48w_rep1","48w_rep2", "48w_rep3")
for(use_plate in all_reps){
  for(use_attr in c("confluence", "GFP.obj.count")){  
    dat <- read.csv(fname_input_data,sep="\t")
    dat <- dat[,str_starts(colnames(dat),"meta") | str_starts(colnames(dat),use_attr)]
    colnames(dat) <- dat[1,]
    dat <- dat[-1,]
    dat <- dat[dat$plate==use_plate,-c(1:2)] #beep! we lost stuff a gene in some analysis
    dat <- melt(dat, id.vars = c(1,2))
    colnames(dat) <- c("time","replicate","gene","value")
    #dat <- dat[!is.na(dat$value),]
    dat$value <- as.double(dat$value)
    dat$time <- as.double(dat$time)
    dat <- dat[dat$value!="NA" & dat$value!="" & !is.na(dat$value),]
    

    ### Normalize by first time point value here?
    first_time_dat <- merge(dat,sqldf::sqldf("select min(time) as time, gene from dat group by gene"))
    normalize_by <- sqldf::sqldf("select gene, avg(value) as first_value from first_time_dat group by gene")
    dat <- merge(dat, normalize_by)
    dat$norm_value <- dat$value/dat$first_value
    
    
    #Compute mean and sd for each time point
    datsum <- as.data.frame(dat %>% group_by(time,gene) %>%  summarize(avg=mean(norm_value), sd=sd(norm_value)))
    datsum$se <- datsum$sd/sqrt(3)

    #GFP count/area makes no sense to show for no virus
    if(use_attr!="confluence"){
      datsum <- datsum[datsum$gene!="no virus",]  
    }
    
    the_xlab <- str_replace(use_plate, "48w_rep","Rep #")
    
    
    p1 <- ggplot(datsum, aes(time, avg, group=gene)) + 
      scale_color_manual(name = "", values = gene_colors) +
      scale_fill_manual(name = "", values = gene_colors) +
      geom_line(mapping = aes(color=gene)) +
      geom_ribbon(aes(y = avg, ymin = avg - se, ymax = avg + se, fill = gene), alpha = .2) +
      theme_bw() + 
      xlab(paste("Time (h)", the_xlab)) + 
      ylab(paste("Normalized", use_attr))  
    allplot[[paste(use_plate, use_attr)]] <- p1
  }
}
ptot <- egg::ggarrange(plots=allplot, nrow=length(all_reps))
ptot
ggsave(plot=ptot, "growth_curves.svg", width = 7, height = 7)





#################
#################  Statistically test the ratio of final cells vs initial cells
#################



all_fits <- NULL
all_attr <- c("confluence", "GFP.obj.count")
for(use_attr in all_attr){ 
  
  ####### Get data in long form for this measurement type
  dat <- read.csv(fname_input_data,sep="\t")
  dat <- dat[,str_starts(colnames(dat),"meta") | str_starts(colnames(dat),use_attr)]
  colnames(dat) <- dat[1,]
  dat <- dat[-1,-2]
  dat <- melt(dat, id.vars = c(1,2,3))
  colnames(dat) <- c("plate","time","replicate","gene","value")
  dat$value <- as.double(dat$value)
  dat$time <- as.double(dat$time)
  dat$attr <- use_attr
  dat$to_fit <- paste(dat$plate, dat$replicate, dat$gene)
  dat <- dat[dat$value!="NA" & dat$value!="" & !is.na(dat$value),]
  
  dat <- dat[!(dat$attr=="GFP.obj.count" & dat$gene=="no virus"),]
  

  ####### Compare each tech rep
  for(for_to_fit in dat$to_fit){
    subdat <- dat[dat$to_fit==for_to_fit,]
    
    #Compare first and last    
    final_value <- subdat$value[subdat$time==max(subdat$time)]
    first_value <- subdat$value[subdat$time==min(subdat$time)]
    the_coef <- log10(final_value/first_value)
    
    ## Add to list of models
    onefit <- data.frame(
      coef=the_coef,
      plate=subdat$plate[1],
      gene=subdat$gene[1],
      attr=subdat$attr[1]
    ) 
    all_fits <- rbind(all_fits, onefit)
  }
}


#Decide what to keep and perform averaging.
all_reps <- c("48w_rep1","48w_rep2", "48w_rep3")
keep_fits <- all_fits[all_fits$plate %in% all_reps,]
allfits_avg <- sqldf::sqldf("select gene, attr, plate, avg(coef) as coef from keep_fits group by gene, attr, plate")
allfits_avg_ntc <- sqldf::sqldf("select attr, plate, coef as ntc_coef from allfits_avg where gene='NTC'")
allfits_avg <- merge(allfits_avg, allfits_avg_ntc)

#Difference in GR
allfits_avg$delta_coef <-  allfits_avg$coef - allfits_avg$ntc_coef


all_gene_stat <- data.frame()
for(use_attr in all_attr){ 
  print(paste("==============",use_attr,"========================="))
  fit_for_attr <- allfits_avg[allfits_avg$attr==use_attr,]
  
  ta <- anova(lm(delta_coef ~ gene, fit_for_attr[fit_for_attr$gene!="no virus",])) #Only compare to gene KOs
  print(paste("anova global P value",ta$`Pr(>F)`[1]))
  
  print("per gene pval, ttest")  
  for(for_gene in unique(fit_for_attr$gene)){
    subdat <- fit_for_attr[fit_for_attr$gene==for_gene,]
    if(use_attr=="GFP.obj.count"){
      subdat <- subdat[subdat$plate!="48w_rep3",]
    }
    tt <- t.test(subdat$coef, subdat$ntc_coef, paired = TRUE)
    print(paste(for_gene,tt$p.value))
    
    all_gene_stat <- rbind(
      all_gene_stat,
      data.frame(
        attr=use_attr,
        gene=for_gene,
        p.value=tt$p.value,
        lower95=tt$conf.int[1],
        upper95=tt$conf.int[2]
      ))
  }
}
all_gene_stat <- all_gene_stat[!is.na(all_gene_stat$p.value),]
all_gene_stat


############ Plotting

use_attr <- "confluence"
fit_for_attr <- allfits_avg[allfits_avg$attr==use_attr,]
p1 <- ggplot(fit_for_attr[fit_for_attr$gene!="NTC",]) + 
  coord_flip() + 
  geom_errorbar(data = all_gene_stat[all_gene_stat$attr==use_attr, ],
                aes(gene,ymin = lower95, ymax = upper95), width=.1, color="gray") +
  geom_point(aes(gene, delta_coef, color=plate)) + 
  xlab("")+
  ylab("Δgrowth vs NTC: Confluence") + 
  theme_bw() +
  theme( 
    panel.grid.major.y = element_blank() ,
    panel.grid.major.x = element_line( size=.3, color="gray" ) 
  )
p1  

use_attr <- "GFP.obj.count"
fit_for_attr <- allfits_avg[allfits_avg$attr==use_attr,]
p2 <- ggplot(fit_for_attr[fit_for_attr$gene!="NTC" & fit_for_attr$plate!="48w_rep3",]) + 
  geom_errorbar(data = all_gene_stat[all_gene_stat$attr==use_attr, ],
                aes(gene,ymin = lower95, ymax = upper95), width=.1, color="gray") +
  geom_point(aes(gene, delta_coef, color=plate)) + 
  coord_flip() + 
  xlab("")+
  ylab("Δgrowth vs NTC: GFP object count")+
  theme_bw() +
  theme( 
    panel.grid.major.y = element_blank() ,
    panel.grid.major.x = element_line( size=.3, color="gray" ) 
  )
p2

ptot <- egg::ggarrange(p1,p2)
ggsave("fitted_growth_rate.svg", ptot, width = 5, height = 3)

all_gene_stat


##################### Plot it for main figure
use_attr <- "confluence"
fit_for_attr <- allfits_avg[allfits_avg$attr==use_attr,]
ptot <- ggplot(fit_for_attr[fit_for_attr$gene!="NTC",]) + 
  coord_flip() + 
  geom_errorbar(data = all_gene_stat[all_gene_stat$attr==use_attr, ],
                aes(gene,ymin = lower95, ymax = upper95), width=.1, color="gray") +
  geom_point(aes(gene, delta_coef, color=plate)) + 
  xlab("")+
  ylab("Δgrowth vs NTC: Confluence") + 
  theme_bw() +
  theme( 
    panel.grid.major.y = element_blank() ,
    panel.grid.major.x = element_line( size=.3, color="gray" ) 
  )
ptot 
ggsave("fitted_growth_rate_main.svg", ptot, width = 4, height = 2)
