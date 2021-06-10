setwd("C:/Users/黄捷/Desktop")
source("D:/scripts/compareP.f.R")

png("compareP.png", w=1200, h=400)
par(mfrow=c(1,3), mai=c(1,1,0.5,0.5))
compareP(
  #f1="D:/files/posCtrls/cad.txt", f1_name="Known", f1_snp="SNP", f1_ea="EA", f1_nea=NA, f1_eaf=NA, f1_beta="BETA", f1_p="P",
  f1="D:/projects/001students/bbc/additive/bb_LDL.gwas.gz", f1_name="UKB", f1_snp="SNP", f1_ea="A1", f1_nea="A2", f1_eaf="A1_FREQ", f1_beta="BETA", f1_p="P",
  f2="D:/data/gwas/bbj/bbj.ldl.gz", f2_name="BBJ", f2_snp="rsids", f2_ea="alt", f2_nea="ref", f2_eaf="maf", f2_beta="beta", f2_p="pval"
)
dev.off()
