setwd("C:/Users/jiehu/Desktop")
source("D:/scripts/library/compareB.f.R")

pdf("compareB.pdf")
par(mfrow=c(2,2), mai=c(1,1,0.5,0.5))
compareP(
  f1="D:/data/gwas/posCtrls/cad.w.t2d.txt", f1_name="CAD", f1_snp="SNP", f1_ea="EA", f1_nea="NEA", f1_eaf="EAF", f1_beta="BETA", f1_se="SE", f1_p="P",
  f2="D:/data/gwas/posCtrls/cad.w.t2d.txt", f2_name="T2D", f2_snp="SNP", f2_ea="EA", f2_nea="NEA", f2_eaf="t2d.EAF", f2_beta="t2d.BETA", f2_se="SE", f2_p="t2d.P"
)
dev.off()
