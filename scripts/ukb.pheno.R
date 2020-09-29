setwd("C:/Users/黄捷/Desktop")
#install.packages(c('data.table', 'VennDiagram', 'dplyr', 'tidyverse')) #("devtools',  "Rcpp", "RcppEigen"))
#install_github("gabraham/plink2R/plink2R") # ukbtools
pacman::p_load(data.table, crosswalkr, dplyr, ggplot2, tidyverse, magrittr, survival, survminer, naniar, VennDiagram, corrplot)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read data 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
typed = read.table(gzfile('D:/projects/001UKB/typed/typed.raw.gz','r'), header=T, as.is=T)
mt = read.table(gzfile('D:/projects/001UKB/typed/mt.raw.gz','r'), header=T, as.is=T)
imp = read.table(gzfile('D:/projects/001UKB/imp/imp.hard.raw.gz','r'), header=T, as.is=T)
hap = read.table('D:/projects/001UKB/hap/apoe.hap.txt', header=T, as.is=T)
wes = read.table('D:/projects/001UKB/wes/apoe.wes.txt', header=T, as.is=T)
names(typed)=gsub("typed.IID", "IID", paste0("typed.", names(typed)))
names(mt)=gsub("mt.IID", "IID", paste0("mt.", names(mt)))
names(imp)=gsub("imp.IID", "IID", paste0("imp.", names(imp)))
names(wes)=gsub("wes.IID", "IID", paste0("wes.", names(wes)))
gen0 = Reduce(function(x,y) merge(x,y,by="IID",all=T), list(typed, mt, imp, hap, wes))
gen0[gen0=="00"] = NA
gen = gen0 %>%
	mutate( 
	abo.O = ifelse( (imp.abo.O1.rs8176719_TC==0 | imp.abo.O2.rs41302905_T==2 | (imp.abo.O1.rs8176719_TC==1 & imp.abo.O2.rs41302905_T>0)), 2, 
		ifelse( (imp.abo.O1.rs8176719_TC==1 | imp.abo.O2.rs41302905_T==1), 1, 0)),
	abo = ifelse(abo.O==2, "OO", 
		ifelse(abo.O==1 & imp.abo.AB.rs8176746_T==0, "OA", 
			ifelse(abo.O==1, "OB", 
				ifelse(imp.abo.AB.rs8176746_T==2, "BB", 
					ifelse(imp.abo.AB.rs8176746_T==1, "AB", "AA"))))),
	sp1 = ifelse(imp.sp1.Z.rs28929474_T==2, "ZZ", 
		ifelse (imp.sp1.S.rs17580_A ==2, "SS", 
			ifelse( (imp.sp1.Z.rs28929474_T==1 & imp.sp1.S.rs17580_A==1), "SZ", 
				ifelse( (imp.sp1.Z.rs28929474_T==1 & imp.sp1.S.rs17580_A==0), "MZ", 
					ifelse( (imp.sp1.Z.rs28929474_T==0 & imp.sp1.S.rs17580_A==1), "MS", 
						ifelse( (imp.sp1.Z.rs28929474_T==0 & imp.sp1.S.rs17580_A==0), "MM", NA))))))
	)

source("D:/projects/001UKB/pheno/vip.r")
pnames <- read.table("D:/files/ukb.vip.fields", header=F)
pnames$V2[duplicated(pnames$V2)] #check duplication
pnames$V1 <- paste0("f.", pnames$V1, ".0.0")
phe <- subset(bd, select=grep("f.eid|\\.0\\.0", names(bd))) %>%
	rename(IID=f.eid) %>%
	renamefrom(pnames, V1, V2, drop_extra=F) %>%
	mutate( race = ifelse(grepl("100",ethnicity),1, ifelse(grepl("200",ethnicity),2, ifelse(grepl("300",ethnicity),3, ifelse(grepl("400",ethnicity),4, ethnicity)))),
	race = factor(race, levels=c(1,2,3,4,5,6), labels=c("White","Mixed","Asian","Black","Chinese","Other")),
	age_cat = cut(age, breaks=seq(35,75,5)),
	sex = factor(sex, levels=c(0,1), labels=c("f","m")),
	edu = factor(edu, levels=1:6, labels=c("college", "A-level", "O-level", "CSE", "NVQ", "other")),
	bmi_cat = cut(bmi, breaks=c(10,18.5,25,30,100), labels=c("lean","healthy","overweight","obese")),
	smoke_status = factor(smoke_status, levels=0:2, labels=c("never","previous","current")), 
	alcohol_status = factor(alcohol_status, levels=0:2, labels=c("never","previous","current"))
	)
phe[phe <0] = NA # BE CAREFUL
icd <- read.table("D:/projects/001UKB/pheno/icd.2cols", header=T, as.is=T) 
ICDnames <- read.table("D:/files/ukb.vip.icd10", header=T)
for (i in 1:nrow(ICDnames)) {
	icd[[paste0("icd_",ICDnames$names[i])]] = ifelse( grepl(ICDnames$codes[i], icd$icd10), 1, ifelse(grepl( substring(ICDnames$codes[i],1,1), icd$icd10), NA,0))
}
sqc <- read.table("D:/projects/001UKB/pheno/raw/ukb_sqc_v2.txt", header=T, as.is=T) %>%
	rename(garray=genotyping.array, in.british=in.white.British.ancestry.subset, in.pca=used.in.pca.calculation, aneuploidy=putative.sex.chromosome.aneuploidy)
sqc <- subset(sqc, select=grepl("IID|garray|in\\.|aneuploidy|kinship|excess|PC", names(sqc)))
rel <- read.table("D:/projects/001UKB/pheno/ukb.related", header=T, as.is=T)
#related <- read.table("D:/projects/001UKB/pheno/raw/ukb1941_rel_s488366.dat", header=T, as.is=T)
#rel <- ukb_gen_samples_to_remove(related, dat0$IID, cutoff = 0.0884)
fatIll <- read.table("D:/projects/001UKB/pheno/fatIll.2cols", header=T, as.is=T)
motIll <- read.table("D:/projects/001UKB/pheno/motIll.2cols", header=T, as.is=T)
phe0 = Reduce(function(x,y) merge(x,y,by="IID",all=T), list(phe, icd, sqc, rel, fatIll, motIll))
dat0 = merge(gen, phe0, by="IID")
dat0 = subset(dat0, IID>0, select=!grepl("IID.", names(dat0)))
saveRDS(dat0, file="D:/projects/001UKB/Rdata/ukb.phe.rds")
dat0$sex12 =ifelse(is.na(dat0$sex), 0, ifelse(dat0$sex=="m", 1, 2))
dat = subset(dat0, select=c("IID","IID","ethnicity","age","sex12","smoke_status","bmi","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")) 
write.table(dat, "ukb.cov", na="NA", append=F, quote=F, col.names=T, row.names=F)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PCA with alcohol and lactase [rs4988235, common T/A is lactase persistent]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat=dat0
eth = c(1001, 1002, 1003, 3001, 3002, 3004, 4001, 4002, 5)
eth_str = c("british", "irish", "white_other", "indian", "pakistani", "asian_other", "carribbean", "african", "chinese")
eth_col = c("gray", "blue", "green", "brown", "brown", "brown", "pink", "black", "yellow")
dat$eth_col = NA; for (i in 1:length(eth)) dat$eth_col[dat$ethnicity==eth[i]] = eth_col[i] 
table(dat$alcohol.ALDH2.rs671_A, dat$alcohol.ADH1B.rs1229984_C)
cor(dat$alcohol.ALDH2.rs671_A, dat$alcohol.ADH1B.rs1229984_C)
table(dat$lactase.MCM6.rs4988235)
dat1 = subset(dat, !is.na(eth_col))
dat2 = subset(dat, alcohol.ALDH2.rs671_A==2)
plot(dat1$PC1, dat1$PC2, col=dat1$eth_col)
points(dat2$PC1, dat2$PC2, cex=1, col="red")
#? how to project 1000G PCs to here?