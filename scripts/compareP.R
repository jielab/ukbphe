setwd("C:/Users/黄捷/Desktop")
# Rscript compareP.R a b c; argv <- commandArgs(TRUE); name1 <- argv[1]

name1="posCtrl"; file1="D:/files/posCtrls/height.EUR3290.txt"; cols1=""
name2="UKB"; file2="D:/projects/001students/001leg/gwas/height.gwas.gz"; cols2="SNP CHR POS EA NEA A1_CT EAF N BETA SE Z P" 

if(grepl(".gz$",file1)) dat1 <- read.table(gzfile(file1,'r'), header=T, as.is=T) else dat1 <- read.table(file1, header=T, as.is=T)
if(grepl(".gz$",file2)) dat2 <- read.table(gzfile(file2,'r'), header=T, as.is=T) else dat2 <- read.table(file2, header=T, as.is=T)
if (cols1 !="") {names(dat1) <- unlist(strsplit(cols1," "))}; dat1$logP <- -log10(dat1$P)
if (cols2 !="") {names(dat2) <- unlist(strsplit(cols2," "))}; dat2$logP <- -log10(dat2$P)
# eval(parse(text=paste0("dat1 <- dat1[, c(", cols1, ")]")))
names(dat1) <- paste0( names(dat1), "_1")
names(dat2) <- paste0( names(dat2), "_2")
dat <- merge(dat1, dat2, by.x="SNP_1", by.y="SNP_2"); dat <- na.omit(dat)
for (var in c('EAF_1', 'EAF_2', 'BETA_1', 'BETA_2', 'P_1', 'P_2')) {
  dat[[var]] <- as.numeric(dat[[var]])
}
#dat <- subset(dat, (EA_1==EA_2 & NEA_1==NEA_2) | (EA_1==NEA_2 & NEA_1==EA_2))
dat[which(dat$EA_1 != dat$EA_2), "EAF_1"]  <- 1- dat[which(dat$EA_1 != dat$EA_2), "EAF_1"]
dat[which(dat$EA_1 != dat$EA_2), "BETA_1"] <- 0- dat[which(dat$EA_1 != dat$EA_2), "BETA_1"]
# dat <- subset(dat, P_1 <0.05 & P_2<0.05); # subset(dat, abs(EAF_1 - EAF_2) <0.4)


pdf(paste0(name1,".",name2,".pdf"), w=10, h=10)
par(mfrow=c(2,2), mai=c(1,1,0.5,0.5))
for (var in c('BETA', 'EAF', 'logP')) {
    X <- dat[[paste0(var,'_1')]]
    Y <- dat[[paste0(var,'_2')]]
    plot(X , Y, main=paste0(var, " (N=", nrow(dat), ")"), xlab=name1, ylab=name2, cex=2, cex.lab=2, cex.axis=2, cex.main=2, font=2, font.lab=2)
}
dev.off()
