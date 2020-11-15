setwd("C:/Users/黄捷/Desktop")
source("D:/scripts/mhplot.f.R")

png("mhplot.png", w=1200, h=1000)
par(mfrow=c(2,1), mai=c(0.6,1,0.6,0)) 
mhplot(trait='CVD', gwas='cvd.gwas.gz', ylim_t=0, poscon_f='cad.txt', dibiao="ukb.chrom.pos.b37")
mhplot(trait='COPD', gwas='copd.gwas.gz', maf='A1_FREQ', ylim_t=3, dibiao="ukb.chrom.pos.b37")
dev.off()
