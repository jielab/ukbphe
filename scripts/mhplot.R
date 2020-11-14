setwd("C:/Users/黄捷/Desktop")
source("D:/scripts/mhplot.f.R")

png("mhplot.png", w=1200, h=1000)
par(mfrow=c(2,1), mai=c(0.6,1,0.6,0)) 
mhplot(trait='cerebral', gwas='D:/projects/001students/001cvd/gwas/cerebro.gwas.gz', ylim_t=3)
mhplot(trait='hypertensive', gwas='D:/projects/001students/001cvd/gwas/hypertensive.gwas.gz', maf='A1_FREQ', poscon_f='D:/files/posCtrls/CAD.txt', ylim_t=3)
dev.off()
