
library(ggplot2)

################################################################################
########## Plot simulation: Variance, function of cells and reads ##############
################################################################################


#p_pickgene=1/20000 # typical p value; one grna per gene

#p_pickgene=1/4000 # typical p value; one grna per gene -- kinome lib

p_pickgene=1/4000/4 # typical p value; one grna per gene -- kinome lib   4x guides/gene


samplefun <- function(p_pickgene, numcell, totread, numsamp=1000){
  p_pickcell <- rbinom(numsamp, numcell, p_pickgene)/numcell
  est_p <- rbinom(numsamp, totread, p_pickcell)/totread
  est_p
  sd(est_p)
}


totest <- merge(
  data.frame(numcell=(1:10)*c(100000)),
  data.frame(totread=(1:10)*c(100000))
)

totest <- merge(
  data.frame(numcell=(1:100)*c(10000)),
  data.frame(totread=(1:100)*c(10000))
)
totest$sd <- NA

for(i in 1:nrow(totest)){
  totest$sd[i] <- samplefun(p_pickgene, totest$numcell[i], totest$totread[i], numsamp=1000)
}

ggplot(totest, aes(numcell, totread, fill=sd/p_pickgene)) + 
  geom_tile() + 
  theme_bw() +
  xlab("Total number of cells") + 
  ylab("Total number of reads")


ggplot(totest, aes(numcell, totread, z=sd)) + #/p_pickgene
  geom_contour() + 
  theme_bw() +
  xlab("Total number of cells") + 
  ylab("Total number of reads")
ggsave("out/supfig7a.svg", width = 4, height = 4)


# two cases: sequencing depth is rate limiting, or the number of cells are




######### pick out a line: variance vs expected variance for poisson


################################################################################
##################### Simulation of CoV, different cell counts #################
################################################################################


### one
onecurve <- data.frame(
  p_pickgene=(1/4000)*(0:40),
  numcell=100000
)
onecurve$sd <- NA
onecurve$totread <- 1000000
onecurve_1 <- onecurve

### two
onecurve <- data.frame(
  p_pickgene=(1/4000)*(0:40),
  numcell=10000
)
onecurve$sd <- NA
onecurve$totread <- 1000000
onecurve_2 <- onecurve


### three -- infinite cells
onecurve <- data.frame(
  p_pickgene=(1/4000)*(0:40),
  numcell=100000000
)
onecurve$sd <- NA
onecurve$totread <- 1000000
onecurve_3 <- onecurve

### all
onecurve <- rbind(onecurve_1,onecurve_2,onecurve_3)

for(i in 1:nrow(onecurve)){
  onecurve$sd[i] <- samplefun(onecurve$p_pickgene[i], onecurve$numcell[i], onecurve$totread[i], numsamp=100000)
}
#ggplot(onecurve, aes(p_pickgene, sd/p_pickgene, color=totread, group=totread)) + geom_line() + ylab("COV(p_pickgene)")

onecurve$cov <- onecurve$sd / onecurve$p_pickgene


onecurve$grp <- paste0("#Cell:",onecurve$numcell," #Reads:",onecurve$totread)
ggplot(onecurve, aes(p_pickgene, cov, color=grp, group=grp)) + 
  geom_line() + 
  xlab("Abundance") + ylab("COV(Abundance)") + 
  theme_bw()
ggsave("out/supfig7b.svg", width = 6, height = 4)



###Poisson:
#mean is lambda; var is lambda. COV is 1/sqrt(lambda)



################################################################################
##################### PCR variance etc ######################################### new eq
################################################################################

K <- 0.95
mu <- K+1
sigma2 <- (1-mu)**2 * (1-K) + (2-mu)**2*K

f_with_startnum <- function(startnum){
  pcr_est <- data.frame(
    n=1:25
  )
  pcr_est$startnum <- startnum 
  
  pcr_est$var_xn <- mu**(pcr_est$n-1) * (1-mu**pcr_est$n) / (1-mu) * sigma2
  pcr_est$mean_xn <- mu**pcr_est$n * pcr_est$startnum  
  pcr_est
}

toplot <- rbind(f_with_startnum(1), f_with_startnum(2), f_with_startnum(3), f_with_startnum(4), f_with_startnum(32), f_with_startnum(1024))
toplot$sd_xn <- sqrt(toplot$var_xn)
toplot$cov_xn <- sqrt(toplot$var_xn)/toplot$mean_xn

toplot$startnum_f <- factor(paste(toplot$startnum), levels=unique(paste(toplot$startnum)))

ggplot(toplot, aes(mean_xn, cov_xn, color=startnum_f, group=startnum_f)) + 
  geom_point() + 
  geom_line() + 
  scale_x_log10() + 
  scale_y_log10() +
  theme_bw() + 
  xlab("E[X_n]") + 
  ylab("CoV[X_n]") +
  scale_color_discrete(name = bquote(X[1]))
ggsave("out/supfig7c.svg", width = 5, height = 4)





################################################################################
##################### Estimation variance due to PCR  ##########################  hard!
################################################################################

K <- 0.95
mu <- K+1
sigma2 <- (1-mu)**2 * (1-K) + (2-mu)**2*K

f_var <- function(){
  pcr_est <- data.frame(
    startnum = 1:100
  )
  pcr_est$n <- 25
  #  pcr_est$startnum <- startnum 
  
  pcr_est$var_xn <- mu**(pcr_est$n-1) * (1-mu**pcr_est$n) / (1-mu) * sigma2
  pcr_est$mean_xn <- mu**pcr_est$n * pcr_est$startnum  
  pcr_est
}


