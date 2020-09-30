setwd("C:/Users/黄捷/Desktop")
pacman::p_load(data.table, dplyr, ggplot2, tidyverse, magrittr, survival, survminer, naniar, VennDiagram, corrplot)


dat0 <- readRDS("D:/projects/001UKB/Rdata/ukb.phe.rds")
dat = dat0 %>% filter( race=="White" & is.na(related)) 
cvd_names = c("artery", "cerebro", "ischaemic", "pulmonary", "vein", "rheumatic", "other", "unspecified")
dat$cvd_cnt = rowSums(dat[, paste0("icd_",cvd_names)], na.rm=T)
for (trait in cvd_names) {
	dat[[trait]] = ifelse(dat$icd_cvd==0, 0, ifelse(dat[[paste0("icd_",trait)]]==1 & dat$cvd_cnt==1, 1, NA))
}

dat1 <- subset(dat, select=paste0("icd_",cvd_names)); dat1[is.na(dat1)] =0; names(dat1) = cvd_names
dat1.mat <- matrix(as.numeric(as.matrix(dat1)), nrow = nrow(dat1))
dat1.sum <- crossprod(dat1.mat)
colsum <- apply(dat1.mat, 1, sum); dat1.mat[colsum !=1, ] =0
rowsum <- apply(dat1.mat, 2, sum); diag(dat1.sum) <- rowsum
rownames(dat1.sum) <- cvd_names; colnames(dat1.sum) <- cvd_names; dat1.sum
