setwd("C:/Users/黄捷/Desktop")
source("D:/scripts/compareP.f.R")

pdf("compareP.pdf", w=12, h=10)
par(mfrow=c(2,2), mai=c(1,1,0.5,0.5))
compareP(
	name1="posCtrl", 
	file1="D:/files/posCtrls/height.EUR3290.txt", 
	cols1=NA, 
	name2="UKB", 
	file2="D:/projects/001students/001leg/gwas/height.gwas.gz", 
	cols2="SNP CHR POS EA NEA A1_CT EAF N BETA SE Z P",
	plots="BETA EAF logP"
	)
dev.off()



