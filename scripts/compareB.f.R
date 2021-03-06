compareB <-function(
	f1, f1_name, f1_snp, f1_ea, f1_nea, f1_eaf, f1_beta, f1_se, f1_p,
	f2, f2_name, f2_snp, f2_ea, f2_nea, f2_eaf, f2_beta, f2_se, f2_p
	){
	
	if(grepl(".gz$",f1)) dat1 <- read.table(gzfile(f1,'r'), header=T, as.is=T) else dat1 <- read.table(f1, header=T, as.is=T)
	if(grepl(".gz$",f2)) dat2 <- read.table(gzfile(f2,'r'), header=T, as.is=T) else dat2 <- read.table(f2, header=T, as.is=T)
	dat1$POS <- NULL; dat2$POS <- NULL; # avoid "POS" is used for "P"
	
	if (is.na(f1_nea)) {f1_nea="NEA"; dat1$NEA=NA}
	if (is.na(f2_nea)) {f2_nea="NEA"; dat2$NEA=NA}
	if (is.na(f1_eaf)) {f1_eaf="EAF"; dat1$EAF=0}
	if (is.na(f2_eaf)) {f2_eaf="EAF"; dat2$EAF=0}
	if (is.null(dat1[f1_beta]) & !is.null(dat1$OR)) dat1[f1_beta]=log(dat1$OR)
	if (is.null(dat2[f2_beta]) & !is.null(dat2$OR)) dat2[f2_beta]=log(dat2$OR)
	if (is.na(f1_se)) {f1_se="SE"; dat1$SE=NA}
	if (is.na(f2_se)) {f2_se="SE"; dat2$SE=NA}
	if (is.null(dat1[f1_p])) dat1[f1_p]=NA; dat1[f1_p] <- -log10(dat1[f1_p])
	if (is.null(dat2[f2_p])) dat2[f2_p]=NA; dat2[f2_p] <- -log10(dat2[f2_p])
	
	dat1 <- subset(dat1, select=c(f1_snp, f1_ea, f1_nea, f1_eaf, f1_beta, f1_se, f1_p))
	dat2 <- subset(dat2, select=c(f2_snp, f2_ea, f2_nea, f2_eaf, f2_beta, f2_se, f2_p))	
	names(dat1) <- paste0( c("SNP", "EA", "NEA", "EAF", "BETA", "SE", "P"), "_1")
	names(dat2) <- paste0( c("SNP", "EA", "NEA", "EAF", "BETA", "SE", "P"), "_2")
	dat <- merge(dat1, dat2, by.x="SNP_1", by.y="SNP_2")
	for (var in c("EA_1", "NEA_1", "EA_2", "NEA_2")) { dat[[var]] <- toupper(dat[[var]]) }
	if(nrow(dat) > 1000) { dat <- subset(dat, P_1>5) }
	# dat <- na.omit(dat)
	#for (var in c('EAF_1', 'EAF_2', 'BETA_1', 'BETA_2', 'P_1', 'P_2')) { dat[[var]] <- as.numeric(dat[[var]]) }
	#dat <- subset(dat, (EA_1==EA_2 & NEA_1==NEA_2) | (EA_1==NEA_2 & NEA_1==EA_2))
	#dat[which(dat$EA_1 != dat$EA_2), "EAF_2"]  <- 1- dat[which(dat$EA_1 != dat$EA_2), "EAF_2"]
	dat$EAF_2 = ifelse(dat$EA_1 == dat$EA_2, dat$EAF_2, 1- dat$EAF_2)
	dat$BETA_2 = ifelse(dat$EA_1 == dat$EA_2, dat$BETA_2, 0- dat$BETA_2)
	dat$BETA.new_1 = ifelse(dat$BETA_1>0, dat$BETA_1, 0- dat$BETA_1)
	dat$BETA.new_2 = ifelse(dat$BETA_1>0, dat$BETA_2, 0- dat$BETA_2)	
	write.table(dat, paste(f1_name,f2_name,"merged.txt",sep="."), na="NA", append=F, quote=F, col.names=T, row.names=F)
	
	for (var in c("EAF", "P", "BETA", "BETA.new")) { 
		dat1 <- subset(dat, select=c(paste0(var,'_1'), paste0(var,'_2'))) # Var 1 on X-axis
		dat1 <- na.omit(dat1)
		if (nrow(dat1)>0) { # bty="n"
			plot(dat1, main=paste0(var, " (N=", nrow(dat), ")"), xlab=f1_name, ylab=f2_name, pch=20, cex=0.8, cex.axis=1.1, cex.lab=1.2, col=colors()[75], font=2, font.lab=2)
		}
	}
}

