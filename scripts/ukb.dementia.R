setwd("C:/Users/黄捷/Desktop")
pacman::p_load(data.table, dplyr, ggplot2, tidyverse, magrittr, survival, survminer, naniar, corrplot)


### prs ###
prs = read.table(gzfile('D:/projects/001UKB/prs/dementia.prs.gz','r'), header=T, as.is=T) 
#prs = subset(prs, select=grepl("IID$|prs", names(prs)))
#prs.tmp =subset(prs, select=grepl("prs.chr", names(prs)))
#corval=cor(prs.tmp);  round(corval,3);  corval[corval==1] = NA; summary(abs(as.numeric(corval)))
#corrplot(cor(prs.tmp), col=terrain.colors(100), cl.pos="b", tl.pos="d", tl.srt=60)
#corrplot.mixed(cor(prs.tmp), number.cex=1, number.digits=2, tl.cex=1, cl.cex=1)
#prs$prs = rowSums(subset(prs, select=grepl("prs.chr", names(prs))), na.rm=T)


### read data ###
dat0 <- readRDS("D:/projects/001UKB/Rdata/ukb.phe.rds")
dat0 <- subset(dat0, select=grep("IID|array|race|related|age|sex|Ill|bmi|attend|birth|death|smoke|alcohol|PC1$|PC2$|apoe|_ad|_dementia|dep|dis|bc_", names(dat0))) 
dat0 <- merge(dat0, prs, by="IID") %>%
	rename(date_dx = date_dementia) %>%
	mutate(
		dx = ifelse( is.na(date_dx), 0,1),
		fatIll=gsub("-17|-27","99",fatIll), motIll=gsub("-17|-27","99",motIll),
		fatDx =ifelse(grepl("10", fatIll), 1, ifelse(grepl("^\\d|,\\d", fatIll), 0, NA)),
		motDx =ifelse(grepl("10", motIll), 1, ifelse(grepl("^\\d|,\\d", motIll), 0, NA)),
		parentDx=ifelse(is.na(fatDx) & is.na(motDx), NA, ifelse(fatDx==1 & motDx==1, "both", ifelse(fatDx==1, "fat", ifelse(motDx==1, "mot", "none")))),
		parentDx=factor(parentDx, levels=c("none", "fat", "mot", "both")),
		apoe = factor(apoe, levels=c("e3e3", "e2e2", "e2e3", "e2e4", "e3e4", "e4e4")),
		e4dose = ifelse(apoe=="e4e4", 2, ifelse(grepl("e4", apoe), 1, 0)),
		e4dosef= factor(e4dose, levels=0:2, labels=c("zero", "one", "two")),
		e4yes = ifelse(grepl("e4",apoe), 1, 0)
	)
table(dat0$fatDx); table(dat0$motDx); table(dat0$fatDx, dat0$motDx)
table(dat0$parentDx, dat0$fatDx); table(dat0$parentDx, dat0$motDx); table(dat0$parentDx, useNA="always")
prop.table(table(dat0$dx, dat0$fatDx), 1); prop.table(table(dat0$dx, dat0$motDx), 1)
prop.table(table(dat0$e4dosef, dat0$dx), 1)


### sanity check ###
table(dat0$smoke_ever, dat0$smoke_status); table(dat0$alcohol_status)
summary(dat0$date_dementia); summary(dat0$date_ad)
ggplot(dat0, aes(fill=apoe, y=bc_LDL, x=race)) + geom_bar(position="dodge", stat="identity")
bp <- boxplot(bc_LDL ~ race*apoe, data=dat0, xlab="", ylab="", main="", las=2, col=rainbow(6), font=2); bp$stats
aggregate(bc_LDL ~ race*apoe, data=dat0, FUN=function(x) {round(c(length(x), mean(x), sd(x), quantile(x,probs=c(0,0.5,1))), 2)} )
dat <- subset(dat0, select=grep("rheu|oest", grep("sex|bc_", names(dat0), value=T), invert=T, value=T)); gg_miss_var(dat, facet=sex)


### data preparation ###
dat <- dat0 %>% 
	filter( race=="White" & age >=60 & !grepl("e1",apoe) ) %>% 
	mutate( 
		# alcohol_freq  # (data-field 1558): 1 -> Daily; 2 -> 3-4 a week; 3 -> 1-2 a week; 4 -> 1-3 a month; 5 -> occasionally; 6 -> Never; -3 -> NA
		# alcohol_amount # (data-field 20403): 1 -> 1/2; 2 -> 3/4; 3 -> 5/6; 4 -> 7/8/9; 5 -> 10 or more
		lastday = fifelse(!is.na(date_dx), date_dx, fifelse(!is.na(death_date), death_date, as.Date("2018-02-28"))),
		followyears = (as.numeric(lastday) - as.numeric(date_attend)) / 365.25
	) %>%
	filter( followyears >0 ) %>% 
	mutate (
		# prs=(prs-mean(prs))/sd(prs),
		prs_qt = cut(prs, breaks=quantile(prs, probs=seq(0,1,0.2), na.rm=T), include.lowest=T, labels=paste0("q",1:5)),
		prs_grp = ifelse(prs_qt=="q1", "low", ifelse(prs_qt=="q5", "high", "middle")),
		prs_grp = factor(prs_grp, levels=c("low","middle","high"))
	)


### sanity check again ###
table(dat$dx, useNA="always"); hist(dat$followyears)
sum(dat$followyears) #1,545,433 person-years reported in paper
prop.table(table(dat$prs_grp, dat$dx), 1) # 1.23% vs. 0.65% reported in paper
prop.table(table(dat$apoe, dat$dx), 1) 
prop.table(table(dat$dx, dat$apoe, dat$smoke_status, dnn=c("AD", "apoe","smoke")), 2)
aggregate(prs ~ apoe, data=dat, mean)
summary( lm( prs ~ apoe, data=dat ))
summary( glm(as.formula(paste0("dx ~ age +sex+PC1+PC2", paste(paste("+", grep("rheu|oest", grep("bc_", names(dat), value=T), invert=T, value=T)), collapse=" "))), data=dat, na.action=na.exclude))
# dementia vs. depression (field 20126)
prop.table(table(dat$sex, dat$icd_depress),1) #0.none, 1:bipolar_1; 2.biplor_2; 3.severe; 4.moderate; 5.rare
prop.table(table(dat$dx, dat$icd_depress),1)
prop.table(table(dat$dx, dat$depress),1)
summary(glm(dx ~ age + sex + depress + icd_depress, data=dat))


### survival analysis ### depress * smoke + alcohol,  
surv.obj <- Surv(time=dat$followyears, event=dat$dx)
#ggsurvplot(survfit(surv.obj ~ apoe, data=dat), ylim=c(0.9,1)) 
summary(cox.fit <- coxph(surv.obj ~ age+sex + smoke_status+alcohol_status + prs+e4yes + parentDx + e4yes * smoke_status, data=dat))
ggforest(cox.fit, data=dat)
summary(cox.fit <- coxph(as.formula(paste0("surv.obj ~ age + sex + smoke_current + alcohol_current + prs_grp + apoe + icd_depress ", paste(paste("+", grep("rheu|oest", grep("bc_", names(dat), value=T), invert=T, value=T)), collapse=" "))), data=dat))
ggforest(cox.fit, data=dat)

