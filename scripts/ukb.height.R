setwd("C:/Users/黄捷/Desktop")
library(data.table); library(dplyr); library(tidyverse); library(magrittr)
library(caret); library(shape); library(glmnet); library(ncvreg); library(corrplot)


### read data ###
dat0 <- readRDS("D:/projects/001UKB/Rdata/ukb.phe.rds")
dat0 <- subset(dat0, select=grepl("IID|race|related|age|sex|attend|birth|death|PC|prs|bmi|height|asthma|copd|icd|cvd", names(dat0))) %>% 
	filter (race=="White" & is.na(related) & height>=140 & height <=210 & height_sitting>=70 & height_sitting <=110 ) %>% 
	mutate(
	leg=height - height_sitting, height_r = height_sitting / leg 
	)
prs <- read.table(gzfile('D:/projects/001UKB/prs/height.2014.prs.gz','r'), header=T, as.is=T)
prs <- subset(prs, select=grepl("IID|PRS", names(prs)))
prs$prs = rowSums(subset(prs, select=grepl("prs.chr", names(prs))), na.rm=T)
prs =subset(prs, select=c("IID", "prs"))
dat <- merge(dat0, prs, by="IID")
# sanity check
m = subset(dat, sex=="m")
m$trait=m$height_sitting
png("jie.png", h=1200, w=800); par(mfrow=c(3,2))
hist(m$trait, breaks=50)
hist(m$trait, breaks=100)
plot(density(m$trait))
par(pty="s"); qqnorm(m$trait); qqline(m$trait)
dev.off()


### height Mendelian Randomization ###
dat$height_pred <- predict.lm(lm(height ~ prs, data=dat), na.action=na.pass)
df <- NULL
for (p in grep("cvd_|icd_", names(dat), value=T)) {
	if ( length(unique(!is.na(dat[[p]]))) > 2 ) {
		tryCatch( lm <- summary(lm(dat[[p]] ~ height_pred +age+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data=dat)), warning=function(w){print(paste("WARN",p));NaN} )
	} else {
		tryCatch( lm <- summary(glm(dat[[p]] ~ height_pred +age+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data=dat, family=binomial)), warning=function(w){print(paste("WARN",p));NaN} )
	}
	df <- rbind(df, c( p, paste0( signif(lm$coef[2,4],2), "|", signif(lm$coef[2,1],2) ) ))
}

