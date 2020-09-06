# ukbphe
UKB Phenotype process
1. ICD 这样的指标，包含了很多不同时间的时间点，量很大，建议分开来处理。下面的命令先生成 icd.tab文件，然后将 icd.tab 文件整合为两列，便于读入R。
	ukbconv ukb42156.enc_ukb r -s42170 -oicd # -s 表示提取一个data-field
	sed -i ‘s/”//g’ icd.tab
	cat icd.tab | sed -e 's/\tNA//g' -e 's/\t/,/2g' | awk '{ if(NR==1) print "IID\ticd10"; else if (NF==1) print $1 "\tNA"; else print $0 }' > icd.2cols


3. 从上述的 ICD.2cols.txt文件里面提取某一个变量，比如 bipolar（对应的ICD-10代码F31），用R读入数据后，用第一行命令生成一个0/1/NA变量。
	phe$icd_bipolar = ifelse(“F31”, phe$icd10), 1, ifelse(“F”, phe$icd10), NA, 0))

4.如果需要批量处理很多ICD，先写一个 VIP.icd.txt 文件，第一列是ICD代码，第二列是相对应的变量的名字，比如I350|I35  stenosis，其中“|”表示“或者”。这个文件第一行写上 codes names，然后用下面的R代码批量执行。
	ICDnames <- read.table("ukb.vip.icd10", header=T)
	for (i in 1:nrow(ICDnames)) { 
phe[[paste0("icd_",ICDnames$names[i])]] = ifelse( grepl(ICDnames$codes[i], phe$icd10), 1, ifelse(grepl( substring(ICDnames$codes[i],1,1), phe$icd10), NA,0)) }

 
#1c. UKB ICD-10 Date数据提取

UKB data showcase 里面的First occurrence表型系列已有各种ICD的Date，例如 http://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=131492。如果想自己动手来弄（不建议），可以用下面的通用代码来生成。

1. 用下面代码，生成一个含有ICD-10和相对应日期的文件。
	echo “41270\n41280” > vip.fields.txt
	ukbconv ukb42156.enc_ukb r -ivip.fields.txt -oicd-date
	sed -i ‘s/”//g’ icd-date.tab

2. 提取单个ICD 的Date, 比如COPD  (代码J440，不是J44) 。
	cnt=`head -1 icd-date.tab | awk '{printf NF}'` # 找出列数
	awk -v cn=$cnt  -v co="J440" '{if (NR==1) print "IID", co; else {c=(cn-1)/2; printf $1;  for (i=2; i<=(c+1); i++) { if ($i==co) printf " "$(i+c) } printf "\n"  }}'   icd-date.tab > icd-date.tmp
	awk ‘NF==2’ icd-date.tmp > icd-date.2cols

3. 对于有一个不同的ICD-Date 的表型，比如 dementia 有5个ICD 代码“F00|F01|F02|F03|G30”，可以按照上述方法分别生成5个文件，比如 icdDate.F00.2cols, icdDate.F01.2cols，等。然后在R里面合并，在 R里面找出最小的日期。  
#2. 跑GWAS

1. 首先，研究人员按照上述步骤提取需要研究的表型数据和相关的covariates，比如 age, sex, PCs。一般来说，quantitative的表型数据要adjust for covariates 和转化成正态分布，这个可以在R里面用下面的命令来实现。对于疾病的binary 表型，只需要指明需要adjust 的covarites，即可。
	trait_res = residuals(lm(trait ~ age+sex+PC1+PC2, na.action=na.exclude)
	trait_inv = qnorm((rank(trait_res,na.last="keep")-0.5) / length(na.omit(trait_res)))

2. 目前，基于UKB的GWAS由专人负责跑，所以研究人员只需要把上述的带有表型数据和covariates的文件发给该专人就可以了。




 
#3. GWAS 结果深度分析 

1. 提取显著信号，添加简单的注释
	plink --annotate MY.gwas.txt NA attrib=snp129.attrib.txt ranges=glist-hg19 --border 10 --pfilter 5e-8 --out MY.gwas.top

2. 用上述方法生成一个比较小的文件，用R里面的qqman package，或者我写的mhplot.R代码，绘制Manhattan plot，也可以绘制QQ plot（不太常用）。

3. 从GWAS catalog (https://www.ebi.ac.uk/gwas) ) 寻找已知信号，通过R 的plot()来比较该 GWAS跟已经发表过的信号的EAF和BETA的一致性。

4. 从千人基因组网站（https://www.internationalgenome.org/data）下载基因数据，作为计算SNP 之间 LD r2的参考。点击该页面 Phase 3 对应的VCF，下载所有以ALL开头的文件。可将下载后的VCF文件的名字改短为ALL.chr*.gz 这样的名字。用下面的命令找出每一条染色体上 independent的top SNPs。
	for chr in {1..22}; do
plink --vcf ALL.chr$chr.gz --clump Height.2018.txt --clump-p1 1e-08 --clump-p2 1e-8 --clump-kb 1000 --clump-r2 0 --out chr$chr
	done  
	如果不考虑 LD, 只考虑距离，可以任意指定一个 LDfile同时设置  --clump-r2=0。
	如果跟已知的GWAS信号文件一并处理，可快速找出本GWAS中跟已知的信号是否强相关。

5. 对于有统计显著性的重点locus，可以ZOOM画图http://locuszoom.org/


 
更多GWAS深度分析 ！

为了熟悉 UKB GWAS数据的格式，可以点击 fastgwa.info 上面的任意一个Phenotype，下载一个或几个完整的GWAS数据，进行测试。
1.	GWAS数据的功能性注释，参照post-GWAS analysis pipeline (github.com/Ensembl/postgap).
2.	多个GWAS 之间的 genetic correlation 分析，请使用LDSC (https://github.com/bulik/ldsc)
3.	多基因风险评分PRS，请使用 LDpred2 https://privefl.github.io/bigsnpr/articles/LDpred2.html
4.	因果分析 Mendelian Randomization，请使用GSMR （https://cnsgenomics.com/software/gsmr/）和MendelianRandomization R package


常用搜索网站
1.	TopMed 浏览器 https://bravo.sph.umich.edu
2.	GnomAD 浏览器https://gnomad.broadinstitute.org
