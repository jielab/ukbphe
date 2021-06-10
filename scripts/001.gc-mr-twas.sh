#!/bin/bash

dir=/mnt/d/data/gwas/bbj #/mnt/d/projects/001students/cvd/gwas
traits="MI STEMI NSTEMI cvd artery cerebro ischaemic pulmonary vein allergy_dx asthma copd"
traits_x="bilirubin bmi.female bmi bmi.male calcium ck creatinine crp dbp hba1c hdl ht ldl menarche menopause monocyte na plt pulsep rbc sbp tg ua wbc"
traits_y="arrhythmia asthma breastcancer cad cataract chf copd lungcancer osteoporosis pad ra stroke t2d"
ldsc_dir=/mnt/d/software_lin/ldsc


## 第一步：LDSC (https://github.com/bulik/ldsc) ##
conda activate ldsc
for trait in $traits_x $traits_y; do
	# zcat $trait.gz | sed 's/\t\t/\tNA\t/g' | gzip -f > $trait.$dat.chk # 确保没有确实数据
	python2 $ldsc_dir/munge_sumstats.py --sumstats $dir/bbj.$trait.gz --chunksize 10000 --snp rsids --a1 alt --a2 ref --frq maf --p pval --N 200000 --signed-sumstats beta,0 --merge-alleles $ldsc_dir/hm3.snplist --out $trait
	# python2 $ldsc_dir/munge_sumstats.py --sumstats $dir/$trait.gwas.gz --chunksize 10000 --snp SNP --a1 A1 --a2 A2 --frq A1_FREQ --p P --N-col N --signed-sumstats Z,0 --merge-alleles $ldsc_dir/hm3.snplist --out $trait
done
for trait in $traits_x $traits_y; do
	echo process $trait
    echo $trait $traits_x $traits_y | sed -e 's/ /.sumstats.gz,/g' -e 's/$/.sumstats.gz/' | xargs -n1 -I % /mnt/d/software_lin/ldsc/ldsc.py --rg % --out $trait.rg --ref-ld-chr $ldsc_dir/eur_w_ld_chr/ --w-ld-chr $ldsc_dir/eur_w_ld_chr/
    awk '$1=="Summary" {printf NR}' $trait.rg.log | xargs -n1 -I % awk -v s=% 'FNR >=s' $trait.rg.log | sed 's/.sumstats.gz//g' > $trait.rg.txt
done


## 第二步：GSMR (https://cnsgenomics.com/software/gcta/#GSMR) ##
for trait in $traits_x $traits_y; do
	echo "SNP A1 A2 freq b se p N" > $trait.gcta.txt
#	zcat $dir/$trait.gwas.gz | awk 'NR>1 {if (NF==18) print $3,$6,$7,$11, $15,$16,$18,$14; else print $3,$6,$7,$9, $11,$12,$14,$10}' >> $trait.gcta.txt
	zcat $dir/bbj.$trait.gz | awk 'NR>1 {if (arr[$5] !="Y") print $5, $4,$3,$10, $8,$9,$7, 200000; arr[$5]="Y"}' | fgrep -v NA  >> $trait.gcta.txt
done
for tx in $traits_x; do
	echo "$tx $tx.gcta.txt" > $tx.exposure
    for ty in $traits_y; do
        if [[ $ty != $tx && ! -f $ty.outcome ]]; then
            echo "$ty $ty.gcta.txt" > $ty.outcome
		fi
		echo -e "
		source('/mnt/d/scripts/library/gsmr_plot.r')
		gsmr_data = read_gsmr_data('$tx.$ty.eff_plot.gz')
		pdf('$tx.$ty.pdf')
		par(mai=c(1,1,0.5,0.5))
		plot_gsmr_effect(gsmr_data, '$tx', '$ty', colors()[75])
		dev.off()
		" > $tx.$ty.plot.R
		echo run gsmr on $tx vs. $ty
		gcta64 --bfile /mnt/d/data/hm3/hm3.b37 --gsmr-file $tx.exposure $ty.outcome --gsmr2-beta --gsmr-direction 2 --diff-freq 0.3 --gwas-thresh 5e-8 --effect-plot --out $tx.$ty
		Rscript $tx.$ty.plot.R
		awk -v t=$tx 'NR==1 || $1==t' $tx.$ty.gsmr > $tx.$ty.way1.gsmr
		awk -v t=$tx 'NR==1 || $2==t' $tx.$ty.gsmr > $tx.$ty.way2.gsmr
	done
done
# summarize all data
fgrep Error *log # make sure no error
awk 'NR==1 || FNR>1' *.way1.gsmr > way1.gsmr
awk 'NR==1 || FNR>1' *.way2.gsmr > way2.gsmr


## FUSION TWAS ##
## 安装 TWAS (http://gusevlab.org/projects/fusion/)
dir_tw=/mnt/d/data/twas_data
dir_gt=$dir_tw/GTEx_v7_multi_tissue
dir_ld=$dir_tw/LDREF
fusion=/mnt/d/software_lin/fusion_twas
for trait in RHR T2D; do
for tissue in `ls -d1 $dir_gt/*/ | sed 's/\/$//' | awk -F '/' '{print $NF}' | awk '{printf " "$1}'`; do
for chr in 7; do
    echo now process trait $trait, tissue $tissue, chr $chr
    Rscript $fusion/FUSION.assoc_test.R --sumstats $dir/summary/$trait.sumstats.gz --chr $chr --out $trait.$tissue.chr$chr.txt --weights $dir_gt/$tissue.P01.pos --weights_dir $dir_gt --ref_ld_chr $dir_ld/1000G.EUR.
done
done
done
