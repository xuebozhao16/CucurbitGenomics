pepo climate associations reveal candidate adaptive variants

# RDA
/data/xuebo/pepo/GEA/01_RDA
## SNP
```shell
```bash
###做RDA的时候是不能有缺失的数据的，因此要先挑出20K的数据，之后把缺失的数据去掉，把基因型改成0,1,2
vcftools --gzvcf /data/xuebo/pepo/pop_genetics/00_207_data/final_SNP/SNP_207.final_chr.vcf.gz --max-missing 1 --maf 0.05 --recode --stdout | bgzip -c -@100 > SNP_nomiss.vcf.gz #2896684
WGS --r 0.007 --model file --type random --file  SNP_nomiss.vcf.gz  --out  SNP_nomiss_20k.vcf & #20065
bgzip -c SNP_nomiss_20k.vcf > SNP_nomiss_20k.vcf.gz
tabix SNP_nomiss_20k.vcf.gz
# 提取表头（样本名）并组合为 header 行
(echo -e "CHROM\tPOS\t$(bcftools query -l SNP_nomiss_20k.vcf.gz | paste -sd '\t')" && \
bcftools query -f '%CHROM\t%POS[\t%GT]\n' SNP_nomiss_20k.vcf.gz) > genotypes_with_header.txt
awk 'BEGIN{OFS="\t"} NR==1{print; next} {
  for (i=1; i<=NF; i++) {
    if ($i == "0/0" || $i == "0|0") $i=0;
    else if ($i == "0/1" || $i == "1/0" || $i == "0|1" || $i == "1|0") $i=1;
    else if ($i == "1/1" || $i == "1|1") $i=2;
    else if ($i == "./." || $i == ".|.") $i=0;
    # 保持其他不变
  }
  print
}' genotypes_with_header.txt | cut -f3- > SNP_nomiss_20k.txt
```

## SV
```shell
```bash
bcftools view -e 'GT==".|."' /data/xuebo/pepo/pop_genetics/00_207_data/final_SV/SV20_207_phasing.vcf > SV_nomiss.vcf #212851
WGS --r 0.1 --model file --type random --file SV_nomiss.vcf  --out SV_nomiss_20k.vcf & #21250
grep -v "#" SV_nomiss_20k.vcf | awk '{for(i=10;i<=NF;i++) count[$i]++} END {for(num in count) print num, count[num]}'
grep -v "##" SV_nomiss_20k.vcf | awk '
{
    # Combine the first and second columns with an underscore
    $1 = $1 "_" $2;
    # Remove columns 2 to 9 by skipping them in the output
    $2 = ""; $3 = ""; $4 = ""; $5 = ""; $6 = ""; $7 = ""; $8 = ""; $9 = "";
    # Print the modified first column and columns from 10 onward
    for (i = 1; i <= NF; i++) {
        if ($i != "") printf "%s\t", $i;
    }
    print ""; # Print a newline after each row
}' | sed -e 's/0|0/0/g' -e 's/1|1/2/g' | sed 's/#CHROM_POS/CHROM_POS/g' | cut -f2- > SV_nomiss_20k.txt
```

# Bayenv
/data/xuebo/pepo/GEA/02_Bayenv
## install
```shell
```bash
##############现在的目的是得到和环境相关的SNP的情况，使用软件包BAYENV
https://bitbucket.org/tguenther/bayenv2_public/src/master/bayenv2_manual.pdf 这是用到的软件和说明书
##软件安装
wget -c https://bitbucket.org/tguenther/bayenv2_public/get/edcea648df0f.zip
unzip edcea648df0f.zip
/data/xuebo/software/tguenther-bayenv2_public-edcea648df0f 把这个加进环境变量里面
PGDSpider2-cli -inputfile ch36_10000.vcf -inputformat VCF -outputfile test.bay -outputformat BAYENV
#还需要一个数据转换的安装包
wget -c https://software.bioinformatics.unibe.ch/pgdspider/PGDSpider_3.0.0.0.zip
unzip PGDSpider_3.0.0.0.zip
java -jar PGDSpider3-cli.jar
```

## SNP
```shell
```bash
/data/xuebo/pepo/GEA/02_Bayenv/SNP
####格式转换，现在是只使用pepo的数据，其他三个亚种的数据都不加上,这样群体结构可以控制一下
cat /data/xuebo/pepo/pop_genetics/group/group_NA.txt /data/xuebo/pepo/pop_genetics/group/group_SA.txt /data/xuebo/pepo/pop_genetics/group/group_EU.txt /data/xuebo/pepo/pop_genetics/group/group_Asia.txt /data/xuebo/pepo/pop_genetics/group/group_AR.txt /data/xuebo/pepo/pop_genetics/group/group_AF.txt > group_only_pepo.txt
vcftools --gzvcf /data/xuebo/pepo/pop_genetics/00_207_data/final_SNP/SNP_207.final_chr.vcf.gz --keep ../group_only_pepo.txt --maf 0.05 --max-missing 1 --recode --stdout | bgzip -c > SNP_onlypepo_maf005.vcf.gz &  #2892464
bcftools annotate --set-id '%CHROM\_%POS' -O z -o SNP_onlypepo_maf005.ID.vcf.gz ../SNP_onlypepo_maf005.vcf

###########################################产生geno文件
WGS --r 0.07 --model file --type random --file  SNP_onlypepo_maf005.vcf.gz  --out  SNP_onlypepo_maf005_200k.vcf & #202386
##准备spid文件
grep -v "C_ozarkana" /data/xuebo/pepo/pop_genetics/group/pepo_189.txt | grep -v "C_texana" | grep -v "Australia" > group_only_pepo2.txt
bcftools query -l SNP_onlypepo_maf005_200k.vcf | paste -sd, -
***************vim SNP.spid
# spid-file generated: Fri Nov 26 09:43:15 CST 2021

# VCF Parser questions
PARSER_FORMAT=VCF

# Only output SNPs with a phred-scaled quality of at least:
VCF_PARSER_QUAL_QUESTION=2
# Select population definition file:
VCF_PARSER_POP_FILE_QUESTION=/data/xuebo/pepo/GEA/02_Bayenv/group_only_pepo2.txt
# What is the ploidy of the data?
VCF_PARSER_PLOIDY_QUESTION=DIPLOID
# Do you want to include a file with population definitions?
VCF_PARSER_POP_QUESTION=TRUE
# Output genotypes as missing if the phred-scale genotype quality is below:
VCF_PARSER_GTQUAL_QUESTION=2
# Do you want to include non-polymorphic SNPs?
VCF_PARSER_MONOMORPHIC_QUESTION=FALSE
# Only output following individuals (ind1, ind2, ind4, ...):
VCF_PARSER_IND_QUESTION=GT004,GT010,GT013,GT014,GT015,GT016,GT017,GT018,GT019,GT020,GT021,GT022,GT023,GT024,GT025,GT026,GT027,GT028,GT029,GT030,GT031,GT032,GT033,GT034,GT035,GT036,GT037,GT038,GT039,GT040,GT041,GT042,GT043,GT044,GT045,GT046,GT047,GT048,GT049,GT050,GT051,GT052,GT053,GT054,GT055,GT056,GT057,GT058,GT059,GT060,GT061,GT062,GT063,GT064,GT065,GT066,GT067,GT068,GT069,GT070,GT071,GT072,GT073,GT074,GT075,GT076,GT077,GT078,GT080,GT081,GT084,GT085,GT086,GT087,GT088,GT089,GT090,GT091,GT092,GT093,GT094,GT095,GT096,GT097,GT098,GT099,GT100,GT101,GT102,GT103,GT110,GT111,GT112,GT113,GT114,GT115,GT116,GT117,GT118,GT119,GT120,GT121,GT122,GT123,GT124,GT125,GT126,GT127,GT128,GT129,GT130,GT131,GT132,GT133,GT134,GT135,GT136,GT137,GT138,GT139,GT140,GT141,GT142,GT143,GT144,GT145,GT146,GT147,GT148,GT149,GT150,GT151,GT152,GT153,GT154,GT155,GT156,GT158,GT159,GT160,GT161,GT164,GT165,GT166,GT167,GT168,GT169,GT170,GT171,GT172,GT173,GT174,GT175,GT176,GT177,GT178,GT179,GT190,GT191,GT192,GT193,GT194,GT195,GT196,GT197,GT198,GT199,GT200,GT201,GT202,GT203,GT204,GT205,GT206,GT207
VCF_PARSER_INDEL_QUESTION=FALSE
# Only input following regions (refSeqName:start:end, multiple regions: whitespace separated):
VCF_PARSER_REGION_QUESTION=chr1:1:26126445 chr2:1:18012633 chr3:1:21830965 chr4:1:14687242 chr5:1:16890741 chr6:1:21953437 chr7:1:16309332 chr8:1:18684281 chr9:1:17695715 chr10:1:19487942 chr11:1:14076622 chr12:1:20760685 chr13:1:18348356 chr14:1:15169833 chr15:1:12448478 chr16:1:22437600 chr17:1:17714425 chr18:1:12288125 chr19:1:16558619 chr20:1:19596639
# Output genotypes as missing if the read depth of a position for the sample is below:
VCF_PARSER_READ_QUESTION=2
# Take most likely genotype if "PL" or "GL" is given in the genotype field?
VCF_PARSER_PL_QUESTION=FALSE
# Do you want to exclude loci with only missing data?
VCF_PARSER_EXC_MISSING_LOCI_QUESTION=TRUE

# BAYENV Writer questions
WRITER_FORMAT=BAYENV

# Save sample file
BAYENV_WRITER_SAMPLE_FILE_QUESTION=FALSE
# Assign half missing genotypes (one allele missing) as complete missing?
BAYENV_WRITER_HALF_MISSING_QUESTION=FALSE
# Do you want to save two additional files with used sample and loci names?
BAYENV_WRITER_WRITE_INFO_FILE_QUESTION=FALSE
# Do you want to save an additional sample file with sample sizes?
BAYENV_WRITER_WRITE_SAMPLE_FILE_QUESTION=FALSE
# Save sample/loci names file
BAYENV_WRITER_INFO_FILE_QUESTION=FALSE

##run PGDSpider3
nohup java -jar -Xmx500g -Xms100g /data/xuebo/software/PGDSpider_3.0.0.0/PGDSpider3-cli.jar  -inputfile SNP_onlypepo_maf005_200k.vcf -inputformat VCF -outputfile SNP_onlypepo_maf005.envgenofile -outputformat BAYENV -spid SNP.spid & #404766行/2 = 202383个位点，三个位点没有算

############################################产生matrix文件
/data/xuebo/pepo/GEA/02_Bayenv/SNP/matrix_file
samplesize.txt
33	24	81	17	13	7
##过滤LD，要保证群体之间没有连锁
bcftools annotate --set-id '%CHROM\_%POS' -O z -o SNP_onlypepo_maf005.ID.vcf.gz ../SNP_onlypepo_maf005.vcf
plink --vcf SNP_onlypepo_maf005.ID.vcf.gz --make-bed --out SNP_onlypepo_maf005 &
plink --bfile SNP_onlypepo_maf005 --indep-pairwise 50 10 0.2 --out SNP_onlypepo_maf005_LD &
sed 's/_/\t/' SNP_onlypepo_maf005_LD.prune.in > SNP_onlypepo_maf005_LD.prune.in.LDlist
bcftools view -R SNP_onlypepo_maf005_LD.prune.in.LDlist SNP_onlypepo_maf005.ID.vcf.gz > SNP_onlypepo_maf005.ID_LD.vcf &
nohup java -jar -Xmx500g -Xms100g /data/xuebo/software/PGDSpider_3.0.0.0/PGDSpider3-cli.jar  -inputfile SNP_onlypepo_maf005.ID_LD.vcf -inputformat VCF -outputfile SNP_onlypepo_LD.envgenofile -outputformat BAYENV -spid ../SNP.spid &
nohup bayenv2 -i SNP_onlypepo_LD.envgenofile -s samplesize.txt -p 6 -k 10000 -r 83556 > SNP_LD.matrix &
tail -n 7 SNP_LD.matrix > matrix_SNP

############################################产生env文件
/Users/xuebozhao/Library/CloudStorage/OneDrive-Personal/FeiLab/pepo/Bayenv/207_bioclim_SNP.csv
/data/xuebo/pepo/GEA/02_Bayenv/SNP/env_file
# 读入环境数据
env <- read.csv("207_bioclim_SNP.csv", row.names = 1)   # 你的左边这个表
group <- read.table("/data/xuebo/pepo/GEA/02_Bayenv/group_only_pepo2.txt", header = FALSE)  # 你的右边这个表
colnames(group) <- c('Sample', 'Group')
rownames(group) <- group$Sample
# 保证顺序一致
env <- env[rownames(group), ]
# 加上群体信息
env$Group <- group$Group
# 只取你要的环境变量列，假设你只要 BIO3 和 BIO4
env_selected <- env[, c('BIO3', 'BIO15', 'Group')]
# 聚合：每个群体取均值
env_grouped <- aggregate(. ~ Group, data = env_selected, FUN = mean)
# 去掉Group列，行名设成Group
rownames(env_grouped) <- env_grouped$Group
env_grouped$Group <- NULL
# 对列做z-score标准化
env_scaled <- scale(env_grouped)
# 最后转置：行是变量，列是群体（Bayenv要求）
env_final <- t(env_scaled)
# 保存文件
write.table(env_final, "env_for_bayenv_SNP.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# 读入环境数据
env <- read.csv("207_bioclim_SV.csv", row.names = 1)   # 你的左边这个表
group <- read.table("/data/xuebo/pepo/GEA/02_Bayenv/group_only_pepo2.txt", header = FALSE)  # 你的右边这个表
colnames(group) <- c('Sample', 'Group')
rownames(group) <- group$Sample
# 保证顺序一致
env <- env[rownames(group), ]
# 加上群体信息
env$Group <- group$Group
# 只取你要的环境变量列，假设你只要 BIO3 和 BIO4
env_selected <- env[, c('BIO4', 'BIO17', 'Group')]
# 聚合：每个群体取均值
env_grouped <- aggregate(. ~ Group, data = env_selected, FUN = mean)
# 去掉Group列，行名设成Group
rownames(env_grouped) <- env_grouped$Group
env_grouped$Group <- NULL
# 对列做z-score标准化
env_scaled <- scale(env_grouped)
# 最后转置：行是变量，列是群体（Bayenv要求）
env_final <- t(env_scaled)
# 保存文件
write.table(env_final, "env_for_bayenv_SV.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

###############################################run
/data/xuebo/pepo/GEA/02_Bayenv/SNP/bayenv2_blocks
#按每 10000 行拆分成文件
split -l 10000 ../SNP_onlypepo_maf005.envgenofile --numeric-suffixes=1 --suffix-length=3 --additional-suffix=.txt block_
####vim calc_bf.sh
#!/bin/bash
#just a small bash script to calculate BFs for all SNPs from SNPFILE
#please copy this script into the same directory as bayenv and execute it there
#please see the Bayenv2 manual for details about usage
#make this script executable (chmod +x calc_bf.sh)
#Usage: ./calc_bf.sh <Name of your SNPSFILE> <Name of your ENVFILE> <Name of your MATFILE> <Nuber of populations> <Number of MCMC iterations> <Number of environmental factors>
SNPFILE=$1
ENVFILE=$2
MATFILE=$3
POPNUM=$4
ITNUM=$5
ENVNUM=$6
OUT=$7
split -a 10 -l 2 $SNPFILE $SNPFILE.snp_batch

for f in $SNPFILE.snp_batch*
do
	bayenv2 -i $f -e $ENVFILE -m $MATFILE -k $ITNUM -r $RANDOM -p $POPNUM -n $ENVNUM -t -o $OUT
done
rm -f $SNPFILE.snp_batch*

#稍微测试一下
bash ../calc_bf.sh block_001.txt /data/xuebo/pepo/GEA/02_Bayenv/SNP/env_file/env_for_bayenv_SNP.txt  /data/xuebo/pepo/GEA/02_Bayenv/SNP/matrix_file/matrix_SNP 6 10000 2 Bayenv_bf.out
#################getbf.sh
bash ../calc_bf.sh block_001.txt /data/xuebo/pepo/GEA/02_Bayenv/SNP/env_file/env_for_bayenv_SNP.txt  /data/xuebo/pepo/GEA/02_Bayenv/SNP/matrix_file/matrix_SNP 6 10000 2 Bayenv_bf.out
for i in {001..041}
do
    mkdir block_${i}
    mv block_${i}.txt ./block_${i}
    cd block_${i}
	bash ../calc_bf.sh block_${i}.txt /data/xuebo/pepo/GEA/02_Bayenv/SNP/env_file/env_for_bayenv_SNP.txt  /data/xuebo/pepo/GEA/02_Bayenv/SNP/matrix_file/matrix_SNP 6 10000 2 block_${i}.out &
    cd ..
done
nohup bash getbf.sh > nohupout 2>& 1 &

#########结果整理
wc -l block_*.out.bf ##有不到5000的
#!/bin/bash
for i in {001..041}
do
    cd block_${i}
    awk -F"snp_batch" '
    function base26_to_num(str,    i, n, c) {
        n = 0
        for (i = 1; i <= length(str); i++) {
            c = substr(str, i, 1)
            n = n * 26 + (ord(c) - ord("a"))
        }
        return n + 1  # 从1开始
    }
    function ord(str) {
        return index("abcdefghijklmnopqrstuvwxyz", str) - 1
    }
    {
        split($2, arr, "\t")
        id = base26_to_num(substr(arr[1], 1, 10))
        print id "\t" arr[2] "\t" arr[3]
    }' block_${i}.out.bf > block_${i}.num.bf
    cd ..
done
#!/bin/bash
for i in {001..041}
do
    echo "awk -v blk=$i '{print ((blk-1) * 5000 + \$1)\"\t\"\$2\"\t\"\$3}' block_${i}/block_${i}.num.bf" >> merge.sh
done
bash merge.sh > SNP_block.num.bf #
########挑出位点来
grep -v "#" ../SNP_onlypepo_maf005_200k.vcf | awk '{print $1"\t"$2}' |grep -v "#" | awk '{print NR"\t"$1"\t"$2}' > SNP.pos.txt
awk 'NR==FNR {pos[$1]=$2"\t"$3; next} ($1 in pos) {print $1 "\t" pos[$1] "\t" $2 "\t" $3}' SNP.pos.txt SNP_block.num.bf > SNP_block.merged.bf
awk '{print $2"\t"$3"\t"$4}' SNP_block.merged.bf > SNP_block.temp.bf
awk '{print $2"\t"$3"\t"$5}' SNP_block.merged.bf > SNP_block.prec.bf
```

## SV
```shell
```bash
/data/xuebo/pepo/GEA/02_Bayenv/SV
####格式转换，现在是只使用pepo的数据，其他三个亚种的数据都不加上,这样群体结构可以控制一下
vcftools --gzvcf /data/xuebo/pepo/pop_genetics/00_207_data/final_SV/SV20_207_unphasing.vcf.gz --keep ../group_only_pepo.txt --maf 0.05 --max-missing 1 --recode --stdout | bgzip -c > SV_onlypepo_maf005.vcf.gz &  #182037
gunzip SV_onlypepo_maf005.vcf.gz
awk 'BEGIN{OFS="\t"} /^#/ {print; next} {$4="A"; $5="G"; print}' SV_onlypepo_maf005.vcf > SV_onlypepo_maf005_2.vcf

###########################################产生geno文件
##准备spid文件
***************vim SV.spid
# spid-file generated: Fri Nov 26 09:43:15 CST 2021

# VCF Parser questions
PARSER_FORMAT=VCF

# Only output SNPs with a phred-scaled quality of at least:
VCF_PARSER_QUAL_QUESTION=2
# Select population definition file:
VCF_PARSER_POP_FILE_QUESTION=/data/xuebo/pepo/GEA/02_Bayenv/group_only_pepo2.txt
# What is the ploidy of the data?
VCF_PARSER_PLOIDY_QUESTION=DIPLOID
# Do you want to include a file with population definitions?
VCF_PARSER_POP_QUESTION=TRUE
# Output genotypes as missing if the phred-scale genotype quality is below:
VCF_PARSER_GTQUAL_QUESTION=2
# Do you want to include non-polymorphic SNPs?
VCF_PARSER_MONOMORPHIC_QUESTION=FALSE
# Only output following individuals (ind1, ind2, ind4, ...):
VCF_PARSER_IND_QUESTION=GT004,GT010,GT013,GT014,GT015,GT016,GT017,GT018,GT019,GT020,GT021,GT022,GT023,GT024,GT025,GT026,GT027,GT028,GT029,GT030,GT031,GT032,GT033,GT034,GT035,GT036,GT037,GT038,GT039,GT040,GT041,GT042,GT043,GT044,GT045,GT046,GT047,GT048,GT049,GT050,GT051,GT052,GT053,GT054,GT055,GT056,GT057,GT058,GT059,GT060,GT061,GT062,GT063,GT064,GT065,GT066,GT067,GT068,GT069,GT070,GT071,GT072,GT073,GT074,GT075,GT076,GT077,GT078,GT080,GT081,GT084,GT085,GT086,GT087,GT088,GT089,GT090,GT091,GT092,GT093,GT094,GT095,GT096,GT097,GT098,GT099,GT100,GT101,GT102,GT103,GT110,GT111,GT112,GT113,GT114,GT115,GT116,GT117,GT118,GT119,GT120,GT121,GT122,GT123,GT124,GT125,GT126,GT127,GT128,GT129,GT130,GT131,GT132,GT133,GT134,GT135,GT136,GT137,GT138,GT139,GT140,GT141,GT142,GT143,GT144,GT145,GT146,GT147,GT148,GT149,GT150,GT151,GT152,GT153,GT154,GT155,GT156,GT158,GT159,GT160,GT161,GT164,GT165,GT166,GT167,GT168,GT169,GT170,GT171,GT172,GT173,GT174,GT175,GT176,GT177,GT178,GT179,GT190,GT191,GT192,GT193,GT194,GT195,GT196,GT197,GT198,GT199,GT200,GT201,GT202,GT203,GT204,GT205,GT206,GT207
VCF_PARSER_INDEL_QUESTION=FALSE
# Only input following regions (refSeqName:start:end, multiple regions: whitespace separated):
VCF_PARSER_REGION_QUESTION=chr1:1:26126445 chr2:1:18012633 chr3:1:21830965 chr4:1:14687242 chr5:1:16890741 chr6:1:21953437 chr7:1:16309332 chr8:1:18684281 chr9:1:17695715 chr10:1:19487942 chr11:1:14076622 chr12:1:20760685 chr13:1:18348356 chr14:1:15169833 chr15:1:12448478 chr16:1:22437600 chr17:1:17714425 chr18:1:12288125 chr19:1:16558619 chr20:1:19596639
# Output genotypes as missing if the read depth of a position for the sample is below:
VCF_PARSER_READ_QUESTION=2
# Take most likely genotype if "PL" or "GL" is given in the genotype field?
VCF_PARSER_PL_QUESTION=FALSE
# Do you want to exclude loci with only missing data?
VCF_PARSER_EXC_MISSING_LOCI_QUESTION=TRUE

# BAYENV Writer questions
WRITER_FORMAT=BAYENV

# Save sample file
BAYENV_WRITER_SAMPLE_FILE_QUESTION=FALSE
# Assign half missing genotypes (one allele missing) as complete missing?
BAYENV_WRITER_HALF_MISSING_QUESTION=FALSE
# Do you want to save two additional files with used sample and loci names?
BAYENV_WRITER_WRITE_INFO_FILE_QUESTION=FALSE
# Do you want to save an additional sample file with sample sizes?
BAYENV_WRITER_WRITE_SAMPLE_FILE_QUESTION=FALSE
# Save sample/loci names file
BAYENV_WRITER_INFO_FILE_QUESTION=FALSE

##run PGDSpider3
nohup java -jar -Xmx500g -Xms100g /data/xuebo/software/PGDSpider_3.0.0.0/PGDSpider3-cli.jar  -inputfile SV_onlypepo_maf005_2.vcf -inputformat VCF -outputfile SV_onlypepo_maf005.envgenofile -outputformat BAYENV -spid SV.spid &

############################################产生matrix文件
/data/xuebo/pepo/GEA/02_Bayenv/SV/matrix_file
cp /data/xuebo/pepo/GEA/02_Bayenv/SNP/matrix_file/matrix_SNP ./matrix_SV

############################################产生env文件
/data/xuebo/pepo/GEA/02_Bayenv/SV/env_file
cp /data/xuebo/pepo/GEA/02_Bayenv/SNP/env_file/env_for_bayenv_SV.txt ./

###############################################run
/data/xuebo/pepo/GEA/02_Bayenv/SV/bayenv2_blocks
#按每 10000 行拆分成文件
split -l 10000 ../SV_onlypepo_maf005.envgenofile --numeric-suffixes=1 --suffix-length=3 --additional-suffix=.txt block_
####vim calc_bf.sh
#!/bin/bash
#just a small bash script to calculate BFs for all SNPs from SNPFILE
#please copy this script into the same directory as bayenv and execute it there
#please see the Bayenv2 manual for details about usage
#make this script executable (chmod +x calc_bf.sh)
#Usage: ./calc_bf.sh <Name of your SNPSFILE> <Name of your ENVFILE> <Name of your MATFILE> <Nuber of populations> <Number of MCMC iterations> <Number of environmental factors>
SNPFILE=$1
ENVFILE=$2
MATFILE=$3
POPNUM=$4
ITNUM=$5
ENVNUM=$6
OUT=$7
split -a 10 -l 2 $SNPFILE $SNPFILE.snp_batch

for f in $SNPFILE.snp_batch*
do
	bayenv2 -i $f -e $ENVFILE -m $MATFILE -k $ITNUM -r $RANDOM -p $POPNUM -n $ENVNUM -t -o $OUT
done
rm -f $SNPFILE.snp_batch*

#稍微测试一下
bash calc_bf.sh block_001.txt /data/xuebo/pepo/GEA/02_Bayenv/SV/env_file/env_for_bayenv_SV.txt  /data/xuebo/pepo/GEA/02_Bayenv/SV/matrix_file/matrix_SV 6 10000 2 Bayenv_bf.out
#################getbf.sh
for i in {001..037}
do
    mkdir block_${i}
    mv block_${i}.txt ./block_${i}
    cd block_${i}
	bash ../calc_bf.sh block_${i}.txt /data/xuebo/pepo/GEA/02_Bayenv/SV/env_file/env_for_bayenv_SV.txt  /data/xuebo/pepo/GEA/02_Bayenv/SV/matrix_file/matrix_SV 6 10000 2 block_${i}.out &
    cd ..
done
nohup bash getbf.sh > nohupout 2>& 1 &

#########结果整理
wc -l ./block_*/block_*.out.bf ##有不到5000的
#!/bin/bash
for i in {001..037}
do
    cd block_${i}
    awk -F"snp_batch" '
    function base26_to_num(str,    i, n, c) {
        n = 0
        for (i = 1; i <= length(str); i++) {
            c = substr(str, i, 1)
            n = n * 26 + (ord(c) - ord("a"))
        }
        return n + 1  # 从1开始
    }
    function ord(str) {
        return index("abcdefghijklmnopqrstuvwxyz", str) - 1
    }
    {
        split($2, arr, "\t")
        id = base26_to_num(substr(arr[1], 1, 10))
        print id "\t" arr[2] "\t" arr[3]
    }' block_${i}.out.bf > block_${i}.num.bf
    cd ..
done
#!/bin/bash
for i in {001..037}
do
    echo "awk -v blk=$i '{print ((blk-1) * 5000 + \$1)\"\t\"\$2\"\t\"\$3}' block_${i}/block_${i}.num.bf" >> merge.sh
done
bash merge.sh > SV_block.num.bf 
########挑出位点来
grep -v "#" ../SV_onlypepo_maf005_2.vcf | awk '{print $1"\t"$2}' |grep -v "#" | awk '{print NR"\t"$1"\t"$2}' > SV.pos.txt
awk 'NR==FNR {pos[$1]=$2"\t"$3; next} ($1 in pos) {print $1 "\t" pos[$1] "\t" $2 "\t" $3}' SV.pos.txt SV_block.num.bf > SV_block.merged.bf
awk '{print $2"\t"$3"\t"$4}' SV_block.merged.bf > SV_block.temp.bf
awk '{print $2"\t"$3"\t"$5}' SV_block.merged.bf > SV_block.prec.bf
```

## top1%和ramdom的位点
/data/xuebo/pepo/GEA/02_Bayenv/site
```shell
```bash
#因为Bayenv的文件出来了，所以有两组的数据，一组是随机的SNP，另一组是adaptive allele的情况
cp /data/xuebo/pepo/GEA/02_Bayenv/SNP/bayenv2_blocks/SNP_block.temp.bf ./
cp /data/xuebo/pepo/GEA/02_Bayenv/SNP/bayenv2_blocks/SNP_block.prec.bf ./
cp /data/xuebo/pepo/GEA/02_Bayenv/SV/bayenv2_blocks/SV_block.temp.bf ./
cp /data/xuebo/pepo/GEA/02_Bayenv/SV/bayenv2_blocks/SV_block.prec.bf ./
#把位点挑出来
sort -k3,3g SNP_block.temp.bf | tail -n 3000 | cut -f1,2 | sort -V -k1,1 -k2,2 > SNP_temp.top1.site.txt
sort -k3,3g SNP_block.prec.bf | tail -n 3000 | cut -f1,2 | sort -V -k1,1 -k2,2 > SNP_prec.top1.site.txt
sort -k3,3g SV_block.temp.bf | tail -n 3000 | cut -f1,2 | sort -V -k1,1 -k2,2 > SV_temp.top1.site.txt
sort -k3,3g SV_block.prec.bf | tail -n 3000 | cut -f1,2 | sort -V -k1,1 -k2,2 > SV_prec.top1.site.txt

sort -k3,3g SNP_block.temp.bf | head -n 50000 | shuf -n 3000 | cut -f1,2 | sort -V -k1,1 -k2,2 > SNP_temp.random.site.txt
sort -k3,3g SNP_block.prec.bf | head -n 50000 | shuf -n 3000 | cut -f1,2 | sort -V -k1,1 -k2,2 > SNP_prec.random.site.txt
sort -k3,3g SV_block.temp.bf | head -n 50000 | shuf -n 3000 | cut -f1,2 | sort -V -k1,1 -k2,2 > SV_temp.random.site.txt
sort -k3,3g SV_block.prec.bf | head -n 50000 | shuf -n 3000 | cut -f1,2 | sort -V -k1,1 -k2,2 > SV_prec.random.site.txt

###############把VCF文件挑出来
#######top1_SNP
cp /data/xuebo/pepo/GEA/02_Bayenv/SNP/SNP_onlypepo_maf005_200k.vcf ./
bgzip -c SNP_onlypepo_maf005_200k.vcf > SNP_onlypepo_maf005_200k.vcf.gz
tabix SNP_onlypepo_maf005_200k.vcf.gz
bcftools view -R ../SNP_temp.top1.site.txt ../SNP_onlypepo_maf005_200k.vcf.gz > SNP_temp.vcf &
bcftools view -R ../SNP_prec.top1.site.txt ../SNP_onlypepo_maf005_200k.vcf.gz > SNP_prec.vcf &
***计算每个分组文件的maf的值
for i in {"NA","SA","EU","Asia","AR","AF"}
do
	vcftools --vcf SNP_temp.vcf --keep /data/xuebo/pepo/pop_genetics/group/group_${i}.txt  --freq --out SNP_temp.${i}.vcf &
    vcftools --vcf SNP_prec.vcf --keep /data/xuebo/pepo/pop_genetics/group/group_${i}.txt  --freq --out SNP_prec.${i}.vcf &
done
***把每个SNP的maf值搞出来
for i in {"NA","SA","EU","Asia","AR","AF"}
do 
    cat  SNP_temp.${i}.vcf.frq | awk 'NR!=1 {print}' | sed "s/:/\t/g" | cut -f 6,8 | awk '{if ($1>=$2) print $2; else print $1}' > pop_temp_${i}.maf
    cat  SNP_prec.${i}.vcf.frq | awk 'NR!=1 {print}' | sed "s/:/\t/g" | cut -f 6,8 | awk '{if ($1>=$2) print $2; else print $1}' > pop_prec_${i}.maf
done
***把SNP的名称搞出来
cat  SNP_temp.NA.vcf.frq | awk 'NR!=1 {print}' | sed "s/:/\t/g" | cut -f 2 | awk '{print "SNP_temp_"$1}' > pop_temp_head.maf
cat  SNP_prec.NA.vcf.frq | awk 'NR!=1 {print}' | sed "s/:/\t/g" | cut -f 2 | awk '{print "SNP_temp_"$1}' > pop_prec_head.maf
*列转置
j=0
for i in {"head","NA","SA","EU","Asia","AR","AF"}
do
    cat pop_temp_${i}.maf  |tr "\n" "\t"|sed -e 's/\t$/\n/' > pop_temp_${j}.maf2
    cat pop_prec_${i}.maf  |tr "\n" "\t"|sed -e 's/\t$/\n/' > pop_prec_${j}.maf2
    j=$((j+1))
done
*合并
for ((i=0;i<=6;i=i+1))do echo pop_temp_$i.maf2 ;done | xargs -i cat {} >> SNP_temp.top1.txt
for ((i=0;i<=6;i=i+1))do echo pop_prec_$i.maf2 ;done | xargs -i cat {} >> SNP_prec.top1.txt

#######ramdom_SNP
cp /data/xuebo/pepo/GEA/02_Bayenv/SNP/SNP_onlypepo_maf005_200k.vcf ./
bgzip -c SNP_onlypepo_maf005_200k.vcf > SNP_onlypepo_maf005_200k.vcf.gz
tabix SNP_onlypepo_maf005_200k.vcf.gz
bcftools view -R ../SNP_temp.random.site.txt ../SNP_onlypepo_maf005_200k.vcf.gz > SNP_temp.vcf &
bcftools view -R ../SNP_prec.random.site.txt ../SNP_onlypepo_maf005_200k.vcf.gz > SNP_prec.vcf &
***计算每个分组文件的maf的值
for i in {"NA","SA","EU","Asia","AR","AF"}
do
	vcftools --vcf SNP_temp.vcf --keep /data/xuebo/pepo/pop_genetics/group/group_${i}.txt  --freq --out SNP_temp.${i}.vcf &
    vcftools --vcf SNP_prec.vcf --keep /data/xuebo/pepo/pop_genetics/group/group_${i}.txt  --freq --out SNP_prec.${i}.vcf &
done
***把每个SNP的maf值搞出来
for i in {"NA","SA","EU","Asia","AR","AF"}
do 
    cat  SNP_temp.${i}.vcf.frq | awk 'NR!=1 {print}' | sed "s/:/\t/g" | cut -f 6,8 | awk '{if ($1>=$2) print $2; else print $1}' > pop_temp_${i}.maf
    cat  SNP_prec.${i}.vcf.frq | awk 'NR!=1 {print}' | sed "s/:/\t/g" | cut -f 6,8 | awk '{if ($1>=$2) print $2; else print $1}' > pop_prec_${i}.maf
done
***把SNP的名称搞出来
cat  SNP_temp.NA.vcf.frq | awk 'NR!=1 {print}' | sed "s/:/\t/g" | cut -f 2 | awk '{print "SNP_temp_"$1}' > pop_temp_head.maf
cat  SNP_prec.NA.vcf.frq | awk 'NR!=1 {print}' | sed "s/:/\t/g" | cut -f 2 | awk '{print "SNP_temp_"$1}' > pop_prec_head.maf
*列转置
j=0
for i in {"head","NA","SA","EU","Asia","AR","AF"}
do
    cat pop_temp_${i}.maf  |tr "\n" "\t"|sed -e 's/\t$/\n/' > pop_temp_${j}.maf2
    cat pop_prec_${i}.maf  |tr "\n" "\t"|sed -e 's/\t$/\n/' > pop_prec_${j}.maf2
    j=$((j+1))
done
*合并
for ((i=0;i<=6;i=i+1))do echo pop_temp_$i.maf2 ;done | xargs -i cat {} >> SNP_temp.random.txt
for ((i=0;i<=6;i=i+1))do echo pop_prec_$i.maf2 ;done | xargs -i cat {} >> SNP_prec.random.txt

#######top1_SV
cp /data/xuebo/pepo/GEA/02_Bayenv/SV/SV_onlypepo_maf005.vcf ./
bgzip -c SV_onlypepo_maf005.vcf > SV_onlypepo_maf005.vcf.gz
tabix SV_onlypepo_maf005.vcf.gz
bcftools view -R ../SV_temp.top1.site.txt SV_onlypepo_maf005.vcf.gz > SV_temp.vcf &
bcftools view -R ../SV_prec.top1.site.txt SV_onlypepo_maf005.vcf.gz > SV_prec.vcf &
***计算每个分组文件的maf的值
for i in {"NA","SA","EU","Asia","AR","AF"}
do
	vcftools --vcf SV_temp.vcf --keep /data/xuebo/pepo/pop_genetics/group/group_${i}.txt  --freq --out SV_temp.${i}.vcf &
    vcftools --vcf SV_prec.vcf --keep /data/xuebo/pepo/pop_genetics/group/group_${i}.txt  --freq --out SV_prec.${i}.vcf &
done
***把每个SV的maf值搞出来 
for i in {"NA","SA","EU","Asia","AR","AF"}
do 
    cat  SV_temp.${i}.vcf.frq | awk 'NR!=1 {print}' | sed "s/:/\t/g" | cut -f 6,8 | awk '{if ($1>=$2) print $2; else print $1}' > pop_temp_${i}.maf
    cat  SV_prec.${i}.vcf.frq | awk 'NR!=1 {print}' | sed "s/:/\t/g" | cut -f 6,8 | awk '{if ($1>=$2) print $2; else print $1}' > pop_prec_${i}.maf
done
***把SV的名称搞出来
cat  SV_temp.NA.vcf.frq | awk 'NR!=1 {print}' | sed "s/:/\t/g" | cut -f 2 | awk '{print "SV_temp_"$1}' > pop_temp_head.maf
cat  SV_prec.NA.vcf.frq | awk 'NR!=1 {print}' | sed "s/:/\t/g" | cut -f 2 | awk '{print "SV_temp_"$1}' > pop_prec_head.maf
*列转置
j=0
for i in {"head","NA","SA","EU","Asia","AR","AF"}
do
    cat pop_temp_${i}.maf  |tr "\n" "\t"|sed -e 's/\t$/\n/' > pop_temp_${j}.maf2
    cat pop_prec_${i}.maf  |tr "\n" "\t"|sed -e 's/\t$/\n/' > pop_prec_${j}.maf2
    j=$((j+1))
done
*合并
for ((i=0;i<=6;i=i+1))do echo pop_temp_$i.maf2 ;done | xargs -i cat {} >> SV_temp.top1.txt
for ((i=0;i<=6;i=i+1))do echo pop_prec_$i.maf2 ;done | xargs -i cat {} >> SV_prec.top1.txt

#######ramdom_SV
cp /data/xuebo/pepo/GEA/02_Bayenv/SV/SV_onlypepo_maf005.vcf ./
bgzip -c SV_onlypepo_maf005.vcf > SV_onlypepo_maf005.vcf.gz
tabix SV_onlypepo_maf005.vcf.gz
bcftools view -R ../SV_temp.random.site.txt SV_onlypepo_maf005.vcf.gz > SV_temp.vcf &
bcftools view -R ../SV_prec.random.site.txt SV_onlypepo_maf005.vcf.gz > SV_prec.vcf &
***计算每个分组文件的maf的值
for i in {"NA","SA","EU","Asia","AR","AF"}
do
	vcftools --vcf SV_temp.vcf --keep /data/xuebo/pepo/pop_genetics/group/group_${i}.txt  --freq --out SV_temp.${i}.vcf &
    vcftools --vcf SV_prec.vcf --keep /data/xuebo/pepo/pop_genetics/group/group_${i}.txt  --freq --out SV_prec.${i}.vcf &
done
***把每个SV的maf值搞出来
for i in {"NA","SA","EU","Asia","AR","AF"}
do 
    cat  SV_temp.${i}.vcf.frq | awk 'NR!=1 {print}' | sed "s/:/\t/g" | cut -f 6,8 | awk '{if ($1>=$2) print $2; else print $1}' > pop_temp_${i}.maf
    cat  SV_prec.${i}.vcf.frq | awk 'NR!=1 {print}' | sed "s/:/\t/g" | cut -f 6,8 | awk '{if ($1>=$2) print $2; else print $1}' > pop_prec_${i}.maf
done
***把SV的名称搞出来
cat  SV_temp.NA.vcf.frq | awk 'NR!=1 {print}' | sed "s/:/\t/g" | cut -f 2 | awk '{print "SV_temp_"$1}' > pop_temp_head.maf
cat  SV_prec.NA.vcf.frq | awk 'NR!=1 {print}' | sed "s/:/\t/g" | cut -f 2 | awk '{print "SV_temp_"$1}' > pop_prec_head.maf
*列转置
j=0
for i in {"head","NA","SA","EU","Asia","AR","AF"}
do
    cat pop_temp_${i}.maf  |tr "\n" "\t"|sed -e 's/\t$/\n/' > pop_temp_${j}.maf2
    cat pop_prec_${i}.maf  |tr "\n" "\t"|sed -e 's/\t$/\n/' > pop_prec_${j}.maf2
    j=$((j+1))
done
*合并
for ((i=0;i<=6;i=i+1))do echo pop_temp_$i.maf2 ;done | xargs -i cat {} >> SV_temp.random.txt
for ((i=0;i<=6;i=i+1))do echo pop_prec_$i.maf2 ;done | xargs -i cat {} >> SV_prec.random.txt
```

## top5%的基因
/data/xuebo/pepo/GEA/02_Bayenv/gene
```shell
```bash
#因为Bayenv的文件出来了，所以有两组的数据，一组是随机的SNP，另一组是adaptive allele的情况
cp /data/xuebo/pepo/GEA/02_Bayenv/SNP/bayenv2_blocks/SNP_block.temp.bf ./
cp /data/xuebo/pepo/GEA/02_Bayenv/SNP/bayenv2_blocks/SNP_block.prec.bf ./
cp /data/xuebo/pepo/GEA/02_Bayenv/SV/bayenv2_blocks/SV_block.temp.bf ./
cp /data/xuebo/pepo/GEA/02_Bayenv/SV/bayenv2_blocks/SV_block.prec.bf ./
#把位点挑出来
sort -k3,3g SNP_block.temp.bf | tail -n 10118 | cut -f1,2 | sort -V -k1,1 -k2,2 | awk '{print $1"\t"$2-1000"\t"$2+1000}' > SNP_temp.top5.bed
sort -k3,3g SNP_block.prec.bf | tail -n 10118 | cut -f1,2 | sort -V -k1,1 -k2,2 | awk '{print $1"\t"$2-1000"\t"$2+1000}' > SNP_prec.top5.bed
sort -k3,3g SV_block.temp.bf | tail -n 10118 | cut -f1,2 | sort -V -k1,1 -k2,2 | awk '{print $1"\t"$2-1000"\t"$2+1000}' > SV_temp.top5.bed
sort -k3,3g SV_block.prec.bf | tail -n 10118 | cut -f1,2 | sort -V -k1,1 -k2,2 | awk '{print $1"\t"$2-1000"\t"$2+1000}' > SV_prec.top5.bed
#######
for i in "SNP_temp" "SNP_prec" "SV_temp" "SV_prec"
do
    bedtools intersect -a ${i}.top5.bed  -b /data/xuebo/pepo/ref_genome/C39/C39.final.gff3 -wa -wb > selection_${i}_gene.txt
    grep "gene" selection_${i}_gene.txt | awk '{
        match($12, /ID=([^;]+)/, arr);  # 使用正则表达式匹配Name=后面的基因号
        print $1"\t"$2"\t"$3"\t"$7"\t"$8"\t"arr[1];  # 输出第1,2,3,8,9列以及基因号
    }' > selection_${i}_gene2.txt
    join -t $'\t' -1 6 -2 1 <(sort -t $'\t' -k6,6 selection_${i}_gene2.txt) <(sort -t $'\t' -k1,1 /data/xuebo/pepo/ref_genome/C39/C39.final.ahrd) | sort -k2,2n | cut -f1,2,5,6,7 | awk '!seen[$0]++' >  selection_${i}_gene_anno.txt
    cut -f1 selection_${i}_gene_anno.txt > selection_${i}_gene_anno2.txt
done
```

## Pairwise genetic distance 和 Climate distance关系
```shell
```bash
/data/xuebo/pepo/GEA/02_Bayenv/IBS_matrix
cp /data/xuebo/pepo/GEA/02_Bayenv/SNP/SNP_onlypepo_maf005_200k.vcf ./
plink --vcf SNP_onlypepo_maf005_200k.vcf --make-bed --out SNP_data
plink --bfile SNP_data --distance square '1-ibs'  --out SNP_IBS

/data/xuebo/pepo/GEA/02_Bayenv/IBS_matrix
cp /data/xuebo/pepo/GEA/02_Bayenv/SV/SV_onlypepo_maf005.vcf ./
plink --vcf SV_onlypepo_maf005.vcf --make-bed --out SV_data
plink --bfile SV_data --distance square '1-ibs'  --out SV_IBS
```

## Example
```shell
```bash
##装软件
conda install bioconda::ldblockshow

/data/xuebo/pepo/GEA/02_Bayenv/Example
##把这四个文件相同的和独自有的基因的算出来
cp /data/xuebo/pepo/GEA/02_Bayenv/gene/*anno2.txt ./
# 排序去重
sort -u selection_SNP_temp_gene_anno2.txt   > snp_temp.txt
sort -u selection_SNP_prec_gene_anno2.txt   > snp_prec.txt
sort -u selection_SV_temp_gene_anno2.txt    > sv_temp.txt
sort -u selection_SV_prec_gene_anno2.txt    > sv_prec.txt
# 四个文件共有基因（交集）
comm -12 <(comm -12 snp_temp.txt snp_prec.txt) <(comm -12 sv_temp.txt sv_prec.txt) > common_all4.txt
# 所有基因全集
cat snp_temp.txt snp_prec.txt sv_temp.txt sv_prec.txt | sort -u > all_genes.txt
# 合并其他三个文件，用于计算差集
cat snp_prec.txt sv_temp.txt sv_prec.txt | sort -u > others_for_snp_temp.txt
cat snp_temp.txt sv_temp.txt sv_prec.txt | sort -u > others_for_snp_prec.txt
cat snp_temp.txt snp_prec.txt sv_prec.txt | sort -u > others_for_sv_temp.txt
cat snp_temp.txt snp_prec.txt sv_temp.txt | sort -u > others_for_sv_prec.txt
# 差集：每个文件独有的基因
comm -23 snp_temp.txt others_for_snp_temp.txt > only_SNP_temp.txt
comm -23 snp_prec.txt others_for_snp_prec.txt > only_SNP_prec.txt
comm -23 sv_temp.txt others_for_sv_temp.txt   > only_SV_temp.txt
comm -23 sv_prec.txt others_for_sv_prec.txt   > only_SV_prec.txt
# 可选：统计数量
echo "四个文件共有基因数: $(wc -l < common_all4.txt)"
echo "only_SNP_temp: $(wc -l < only_SNP_temp.txt)"
echo "only_SNP_prec: $(wc -l < only_SNP_prec.txt)"
echo "only_SV_temp:  $(wc -l < only_SV_temp.txt)"
echo "only_SV_prec:  $(wc -l < only_SV_prec.txt)"

#############################把这些基因做个注释
###都有的基因
join -t $'\t' -1 1 -2 1 <(sort -t $'\t' -k1,1 common_all4.txt) <(sort -t $'\t' -k1,1 /data/xuebo/pepo/ref_genome/C39/C39.final.ahrd) > common_all4.anno.txt
C39C01G006100（Glutathione S-transferase）
C39C08G007240（Potassium channel SKOR）(这个更加合理：SKOR（Stelar K+ 外向整流器）是一种在植物细胞中发现的钾通道蛋白，特别是在拟南芥中SKOR是电压门控钾通道 Shaker 超家族的成员，在钾离子 (K+) 从根部到植物其他部位的长距离运输中起着至关重要的作用。 具体来说，SKOR 促进 K+ 从星状细胞（木质部周围的细胞）移动到木质部，而木质部是负责在整个植物体内运输水分和营养物质的维管组织。)
C39C15G002610（Peroxidase）
C39C15G003350（ABC transporter）
C39C07G008940（SET domain protein）

###SV temp特有的基因
join -t $'\t' -1 1 -2 1 <(sort -t $'\t' -k1,1 only_SV_temp.txt) <(sort -t $'\t' -k1,1 /data/xuebo/pepo/ref_genome/C39/C39.final.ahrd) > only_SV_temp.anno.txt ##这个太多了，换个更显著点的
#
sort -k3,3g /data/xuebo/pepo/GEA/02_Bayenv/gene/SV_block.temp.bf | tail -n 500 | cut -f1,2 | sort -V -k1,1 -k2,2 | awk '{print $1"\t"$2-1000"\t"$2+1000}' > SV_temp.top05.bed
bedtools intersect -a SV_temp.top05.bed  -b /data/xuebo/pepo/ref_genome/C39/C39.final.gff3 -wa -wb > selection_SV_temp_gene.txt
grep "gene" selection_SV_temp_gene.txt | awk '{
    match($12, /ID=([^;]+)/, arr);  # 使用正则表达式匹配Name=后面的基因号
    print $1"\t"$2"\t"$3"\t"$7"\t"$8"\t"arr[1];  # 输出第1,2,3,8,9列以及基因号
}' > selection_SV_temp_gene2.txt
join -t $'\t' -1 6 -2 1 <(sort -t $'\t' -k6,6  selection_SV_temp_gene2.txt) <(sort -t $'\t' -k1,1 /data/xuebo/pepo/ref_genome/C39/C39.final.ahrd) | sort -k2,2n | cut -f1,2,5,6,7 | awk '!seen[$0]++' | cut -f1 > xxx
##
comm -12 only_SV_temp.txt xxx > xx
join -t $'\t' -1 1 -2 1 <(sort -t $'\t' -k1,1 xx) <(sort -t $'\t' -k1,1 /data/xuebo/pepo/ref_genome/C39/C39.final.ahrd) > only_SV_temp.anno.txt
| 优先级 | 基因 ID                                             | 理由                  |
| --- | ------------------------------------------------- | ------------------- |
| ⭐⭐⭐ | `C39C19G006410`                                   | 冷胁迫响应蛋白（COR family） | cold-regulated 413 plasma membrane protein 2
| ⭐⭐  | `C39C11G011120`                                   | 离子转运，温度胁迫中电解质平衡     | Cation-transporting ATPase
| ⭐⭐  | `C39C12G002660`                                   | 钙信号通路               | calcium-binding protein 39
| ⭐⭐  | `C39C16G001910`                                   | MAPK 信号通路温度感知调控     | mitogen-activated protein kinase kinase kinase 17-like
| ⭐⭐  | `C39C12G007410`                                   | 热响应型 bHLH 转录因子      | transcription factor IBH1
| ⭐⭐  | `C39C14G007110`                                   | NAC 家族 TF，逆境诱导表达    | NAC domain-containing protein 35
| ⭐⭐  | C39C06G011980 – Abscisic stress-ripening protein 2
| ⭐⭐  | C39C08G001670	bHLH147-like TF	转录因子，逆境响应重要角色
| ⭐⭐  | C39C05G012250	VRN1-like TF	与温度感应（春化）相关


###########################################画图 C39C19G006410
grep "C39C19G006410" /data/xuebo/pepo/ref_genome/C39/C39.final.gff3
chr19	EVM	gene	13247982	13249596	.	-	.	ID=C39C19G006410
chr19	EVM	mRNA	13247982	13249596	.	-	.	ID=C39C19G006410.1;Parent=C39C19G006410
chr19	EVM	exon	13249355	13249596	.	-	.	Parent=C39C19G006410.1
chr19	EVM	CDS	13249355	13249596	.	-	0	Parent=C39C19G006410.1
chr19	EVM	exon	13249133	13249240	.	-	.	Parent=C39C19G006410.1
chr19	EVM	CDS	13249133	13249240	.	-	1	Parent=C39C19G006410.1
chr19	EVM	exon	13248970	13249040	.	-	.	Parent=C39C19G006410.1
chr19	EVM	CDS	13248970	13249040	.	-	1	Parent=C39C19G006410.1
chr19	EVM	exon	13247982	13248235	.	-	.	Parent=C39C19G006410.1
chr19	EVM	CDS	13247982	13248235	.	-	2	Parent=C39C19G006410.1
##这个GF的值是3.3202e+00
chr19	13248533	3.3202e+00
bcftools view -r chr19:13242982-13254596 -Oz /data/xuebo/pepo/GEA/02_Bayenv/SV/SV_onlypepo_maf005.vcf.gz -o SV_chr19.vcf.gz
gunzip SV_chr19.vcf.gz
export PERL5LIB=/data/xuebo/miniconda3/lib/site_perl/5.26.2
LDBlockShow -InVCF /data/xuebo/pepo/GEA/02_Bayenv/SV/SV_onlypepo_maf005_2.vcf  -OutPut re5 -Region chr19:1200000:1400000 -OutPng -SeleVar 1
grep -v "##" SV_chr19.vcf | awk 'BEGIN {OFS="\t"} { for (i=10; i<=NF; i++) { if ($i == "0/0") $i = "0"; else if ($i == "1/1") $i = "2"; else if ($i == "0/1") $i = "1"; else if ($i == "1/0") $i = "1" } print }' | cut --complement -f3,6-9 | sed "s/#CHROM/CHROM/g" > SV_chr19.txt ##这个位置在内含子，影响不大，换一个

########################################################C39C12G007410
grep "C39C12G007410" selection_SV_temp_gene2.txt ##chr12	15754030	15756030	15755613	15756056	C39C12G007410
grep "chr12" SV_block.temp.bf | grep "15755030" ##chr12	15755030	2.0186e+00 
grep "chr12" filter_chr12.SV.eff.vcf | grep "15755030" #影响基因下游

########################################################C39C06G011980
grep "C39C06G011980" selection_SV_temp_gene2.txt ##chr6	15577489	15579489	15578294	15578740	C39C06G011980
grep "chr6" SV_block.temp.bf | grep "15578489" ##chr6	15578489	2.2090e+00 ##chr6	15578489	2.2002e-02
grep "chr6" filter_chr6.SV.eff.vcf | grep "15578489" #splice_region_variant&intron_variant|LOW

########################################################C39C08G001670
grep "C39C08G001670" selection_SV_temp_gene2.txt ##chr8	1699393	1701393	1701135	1701779	C39C08G001670
grep "chr8" SV_block.temp.bf | grep "1700393" ##chr8	1700393	1.0193e-01 ##chr8	1700393	2.2216e-02 ##chr8	1700393	2.7387e+00
grep "chr8" filter_chr8.SV.eff.vcf | grep "1700393" #intergenic_region

########################################################C39C05G012250
grep "C39C05G012250" selection_SV_temp_gene2.txt ##chr5	12426396	12428396	12427751	12429243	C39C05G012250
grep "chr5" SV_block.temp.bf | grep "12427396" ##chr5	12427396	2.0554e+00
grep "chr5" filter_chr5.SV.eff.vcf | grep "12427396" #影响基因下游
```

### C39C08G007240
```shell
```bash
###########################################画图 C39C08G007240 (这个基因是大家都有的)
grep "C39C08G007240" /data/xuebo/pepo/ref_genome/C39/C39.final.gff3
chr8	EVM	gene	5231090	5256248	.	-	.	ID=C39C08G007240
chr8	EVM	mRNA	5231090	5256248	.	-	.	ID=C39C08G007240.1;Parent=C39C08G007240
chr8	EVM	exon	5256016	5256248	.	-	.	Parent=C39C08G007240.1
chr8	EVM	CDS	5256016	5256248	.	-	0	Parent=C39C08G007240.1
chr8	EVM	exon	5255257	5255481	.	-	.	Parent=C39C08G007240.1
chr8	EVM	CDS	5255257	5255481	.	-	1	Parent=C39C08G007240.1
chr8	EVM	exon	5255116	5255179	.	-	.	Parent=C39C08G007240.1
chr8	EVM	CDS	5255116	5255179	.	-	1	Parent=C39C08G007240.1
chr8	EVM	exon	5254636	5254975	.	-	.	Parent=C39C08G007240.1
chr8	EVM	CDS	5254636	5254975	.	-	0	Parent=C39C08G007240.1
chr8	EVM	exon	5252592	5252884	.	-	.	Parent=C39C08G007240.1
chr8	EVM	CDS	5252592	5252884	.	-	2	Parent=C39C08G007240.1
chr8	EVM	exon	5247651	5247731	.	-	.	Parent=C39C08G007240.1
chr8	EVM	CDS	5247651	5247731	.	-	0	Parent=C39C08G007240.1
chr8	EVM	exon	5246211	5246309	.	-	.	Parent=C39C08G007240.1
chr8	EVM	CDS	5246211	5246309	.	-	0	Parent=C39C08G007240.1
chr8	EVM	exon	5245445	5245666	.	-	.	Parent=C39C08G007240.1
chr8	EVM	CDS	5245445	5245666	.	-	0	Parent=C39C08G007240.1
chr8	EVM	exon	5241983	5242171	.	-	.	Parent=C39C08G007240.1
chr8	EVM	CDS	5241983	5242171	.	-	0	Parent=C39C08G007240.1
chr8	EVM	exon	5232217	5232292	.	-	.	Parent=C39C08G007240.1
chr8	EVM	CDS	5232217	5232292	.	-	0	Parent=C39C08G007240.1
chr8	EVM	exon	5231530	5231824	.	-	.	Parent=C39C08G007240.1
chr8	EVM	CDS	5231530	5231824	.	-	2	Parent=C39C08G007240.1
chr8	EVM	exon	5231090	5231225	.	-	.	Parent=C39C08G007240.1
chr8	EVM	CDS	5231090	5231225	.	-	1	Parent=C39C08G007240.1
####
selection_SNP_prec_gene2.txt:chr8	5231938	5233938	5231090	5256248	C39C08G007240
selection_SNP_prec_gene2.txt:chr8	5246890	5248890	5231090	5256248	C39C08G007240 
selection_SNP_prec_gene2.txt:chr8	5255657	5257657	5231090	5256248	C39C08G007240
selection_SNP_temp_gene2.txt:chr8	5238907	5240907	5231090	5256248	C39C08G007240
selection_SNP_temp_gene2.txt:chr8	5240723	5242723	5231090	5256248	C39C08G007240
selection_SNP_temp_gene2.txt:chr8	5255657	5257657	5231090	5256248	C39C08G007240
selection_SV_prec_gene2.txt:chr8	5232324	5234324	5231090	5256248	C39C08G007240
selection_SV_prec_gene2.txt:chr8	5247597	5249597	5231090	5256248	C39C08G007240
selection_SV_prec_gene2.txt:chr8	5253544	5255544	5231090	5256248	C39C08G007240
selection_SV_temp_gene2.txt:chr8	5240355	5242355	5231090	5256248	C39C08G007240
##
grep "chr8" SNP_block.prec.bf | grep "5232938" ##chr8	5232938	2.6165e-01
grep "chr8" SNP_block.prec.bf | grep "5247890" ##chr8	5247890	2.9250e-01
grep "chr8" SNP_block.prec.bf | grep "5256657" ##chr8	5256657	1.7817e-01
grep "chr8" SNP_block.temp.bf | grep "5239907" ## chr8	5239907	1.2314e+01
grep "chr8" SNP_block.temp.bf | grep "5241723" ## chr8	5241723	2.1305e-01
grep "chr8" SNP_block.temp.bf | grep "5256657" ## chr8	5256657	2.4021e-01
grep "chr8" SV_block.prec.bf | grep "5233324" ## chr8	5233324	4.8731e-01 ##chr8	5233324	1.8423e-02
grep "chr8" SV_block.prec.bf | grep "5248597" ##chr8	5248597	1.6264e-01
grep "chr8" SV_block.prec.bf | grep "5254544" ##chr8	5254544	3.9401e-01
grep "chr8" SV_block.temp.bf | grep "5241355" ##chr8	5241355	2.3568e-01
##
grep "chr8" filter_chr8.SNP.eff.vcf | grep "5232938" 
grep "chr8" filter_chr8.SNP.eff.vcf | grep "5247890" 
grep "chr8" filter_chr8.SNP.eff.vcf | grep "5256657" 
grep "chr8" filter_chr8.SNP.eff.vcf | grep "5239907" 
grep "chr8" filter_chr8.SNP.eff.vcf | grep "5241723" 
grep "chr8" filter_chr8.SNP.eff.vcf | grep "5256657" #upstream_gene_variant
grep "chr8" filter_chr8.SV.eff.vcf | grep "5233324" 
grep "chr8" filter_chr8.SV.eff.vcf | grep "5248597" 
grep "chr8" filter_chr8.SV.eff.vcf | grep "5254544" 
grep "chr8" filter_chr8.SV.eff.vcf | grep "5241355" 
###
#把SNP和SV的位置合起来
awk 'BEGIN{OFS="\t"} 
{
  if ($0 ~ /^#/) {
    print $0;
  } else {
    for (i=10; i<=NF; i++) {
      split($i, a, ":");
      $i = a[1];
    }
    print $0;
  }
}' /data/xuebo/pepo/GEA/02_Bayenv/SNP/SNP_onlypepo_maf005_200k.vcf > SNP.gt_only.vcf
vcf-concat SNP.gt_only.vcf /data/xuebo/pepo/GEA/02_Bayenv/SV/SV_onlypepo_maf005.vcf > xxx.vcf
bcftools sort xxx.vcf -o SNP_SV.vcf
awk 'BEGIN{OFS="\t"} 
{
  if ($0 ~ /^#/) {
    print $0;
  } else {
    for (i=10; i<=NF; i++) {
      split($i, a, ":");
      $i = a[1];
    }
    print $0;
  }
}' SNP_SV.vcf > SNP_SV_2.vcf
bgzip -c SNP_SV_2.vcf > SNP_SV_2.vcf.gz
tabix SNP_SV_2.vcf.gz
####
bcftools view -r chr8:5229090-5258248 -Oz SNP_SV_2.vcf.gz -o SV_chr8.vcf.gz
gunzip SV_chr8.vcf.gz
awk 'BEGIN{OFS="\t"}{if($0~"^#"){print}else{for(i=10;i<=NF;i++){split($i,a,":");$i=a[1]}print}}' SV_chr8.vcf > xx
export PERL5LIB=/data/xuebo/miniconda3/lib/site_perl/5.26.2
LDBlockShow -InVCF SNP_SV_2.vcf -OutPut re5 -Region chr8:5100000:5400000 -OutPng -SeleVar 1
grep -v "##" xx | awk 'BEGIN {OFS="\t"} { for (i=10; i<=NF; i++) { if ($i == "0/0") $i = "0"; else if ($i == "1/1") $i = "2"; else if ($i == "0/1") $i = "1"; else if ($i == "1/0") $i = "1" } print }' | cut --complement -f3,6-9 | sed "s/#CHROM/CHROM/g" > SV_chr8.txt
```

### C39C14G007110
```shell
```bash
########################################################C39C14G007110
grep "C39C14G007110" selection_SV_temp_gene2.txt ##chr14	6567869	6569869	6567126	6568618	C39C14G007110
grep "chr14" SV_block.temp.bf | grep "6568869" ##chr14	6568869	2.6152e-02 ##chr14	6568869	2.6431e+0
grep "chr14" filter_chr14.SV.eff.vcf | grep "6568869" #影响基因上游,插入长度约68bp
grep "C39C14G007110" /data/xuebo/pepo/ref_genome/C39/C39.final.gff3
chr14	EVM	gene	6567126	6568618	.	-	.	ID=C39C14G007110
chr14	EVM	mRNA	6567126	6568618	.	-	.	ID=C39C14G007110.1;Parent=C39C14G007110
chr14	EVM	exon	6568354	6568618	.	-	.	Parent=C39C14G007110.1
chr14	EVM	CDS	6568354	6568618	.	-	0	Parent=C39C14G007110.1
chr14	EVM	exon	6567945	6568213	.	-	.	Parent=C39C14G007110.1
chr14	EVM	CDS	6567945	6568213	.	-	2	Parent=C39C14G007110.1
chr14	EVM	exon	6567126	6567821	.	-	.	Parent=C39C14G007110.1
chr14	EVM	CDS	6567126	6567821	.	-	0	Parent=C39C14G007110.1
bcftools view -r chr14:6562126-6573618 -Oz /data/xuebo/pepo/GEA/02_Bayenv/SV/SV_onlypepo_maf005.vcf.gz -o SV_chr14.vcf.gz
gunzip SV_chr14.vcf.gz
export PERL5LIB=/data/xuebo/miniconda3/lib/site_perl/5.26.2
LDBlockShow -InVCF /data/xuebo/pepo/GEA/02_Bayenv/SV/SV_onlypepo_maf005_2.vcf  -OutPut re5 -Region chr14:6400000:6700000 -OutPng -SeleVar 1
grep -v "##" SV_chr14.vcf | awk 'BEGIN {OFS="\t"} { for (i=10; i<=NF; i++) { if ($i == "0/0") $i = "0"; else if ($i == "1/1") $i = "2"; else if ($i == "0/1") $i = "1"; else if ($i == "1/0") $i = "1" } print }' | cut --complement -f3,6-9 | sed "s/#CHROM/CHROM/g" > SV_chr14.txt
```


# Genetic load
### SNP_fold4
```shell
```bash
##/data/xuebo/pepo/pop_genetics/01_tree/outgroup 使用的是这个
cp merged_shared_4D_site_addout.vcf.gz /data/xuebo/pepo/load/SNP_fold4
bcftools annotate -x ^FORMAT/GT -O z -o SNP_site_addout.vcf.gz merged_shared_4D_site_addout.vcf.gz
zcat SNP_site_addout.vcf.gz | grep -v "#" | awk '{print NF; exit}' #有219列,watermelon和Cargyrosperma分别在217和218列
zcat SNP_site_addout.vcf.gz | grep -v "#" | awk 'BEGIN {OFS=FS="\t"} {
    if ($217 == "1/1" && $218 == "1/1") {
        $3 = "D"
        print
    } else if ($217 == "0/0" && $218 == "0/0") {
        $3 = "A"
        print
    }
}' > SNP_addout #13274
awk 'BEGIN {OFS=FS="\t"} { NF = NF - 3; print}' SNP_addout > snp
zcat SNP_site_addout.vcf.gz | grep "#" > header.txt #手动把三个outgroup去掉
cat header.txt snp > snp_add_anc.vcf 
##把有缺失的位点去掉
bcftools view -e 'GT=="./."' snp_add_anc.vcf > snp_add_anc_nomiss.vcf #13088
grep -v "##" snp_add_anc.vcf | awk '{print NF; exit}' #216
grep -v "#" snp_add_anc_nomiss.vcf | awk '{for(i=10;i<=NF;i++) count[$i]++} END {for(num in count) print num, count[num]}'
grep -v "#" snp_add_anc.vcf | awk 'BEGIN {A_count=0; D_count=0} {if ($3 == "A") {A_count++} else if ($3 == "D") {D_count++}} END {print "A_count:", A_count, "D_count:", D_count}'  #A_count: 1638 D_count: 11636
grep -v "#" snp_add_anc_nomiss.vcf | awk 'BEGIN {A_count=0; D_count=0} {if ($3 == "A") {A_count++} else if ($3 == "D") {D_count++}} END {print "A_count:", A_count, "D_count:", D_count}'  #A_count: 1580 D_count: 11448
####derive allele count
bcftools query -l snp_add_anc_nomiss.vcf > xx
    grep -v "#" snp_add_anc_nomiss.vcf | \
    awk 'BEGIN {OFS="\t"} {
        if ($3 == "A") {
            for (i=10; i<=NF; i++) {
                if ($i == "1/1") $i = 2;
                else if ($i == "0/1" || $i == "1/0") $i = 1;
                else if ($i == "0/0") $i = 0;
            }
        } else if ($3 == "D") {
            for (i=10; i<=NF; i++) {
                if ($i == "1/1") $i = 0;
                else if ($i == "0/1" || $i == "1/0") $i = 1;
                else if ($i == "0/0") $i = 2;
            }
        }
        for (i=10; i<=NF; i++) {
            printf "%s%s", $i, (i==NF ? "\n" : OFS)
        }
    }' | awk '{ for(i=1; i<=NF; i++) sum[i] += $i } END { for(i=1; i<=NF; i++) print "sample"i"\t", sum[i] }' | paste xx - | sort -k1,1n > fold4.count.txt
```

### SNP_fold0
```shell
```bash
##/data/xuebo/pepo/pop_genetics/01_tree/outgroup2 使用的是这个
cp merged_shared_4D_site_addout.vcf.gz /data/xuebo/pepo/load/SNP_fold0
bcftools annotate -x ^FORMAT/GT -O z -o SNP_site_addout.vcf.gz merged_shared_4D_site_addout.vcf.gz
zcat SNP_site_addout.vcf.gz | grep -v "#" | awk '{print NF; exit}' #有219列,watermelon和Cargyrosperma分别在217和218列
zcat SNP_site_addout.vcf.gz | grep -v "#" | awk 'BEGIN {OFS=FS="\t"} {
    if ($217 == "1/1" && $218 == "1/1") {
        $3 = "D"
        print
    } else if ($217 == "0/0" && $218 == "0/0") {
        $3 = "A"
        print
    }
}' > SNP_addout #16577
awk 'BEGIN {OFS=FS="\t"} { NF = NF - 3; print}' SNP_addout > snp
zcat SNP_site_addout.vcf.gz | grep "#" > header.txt #手动把三个outgroup去掉
cat header.txt snp > snp_add_anc.vcf 
##把有缺失的位点去掉
bcftools view -e 'GT=="./."' snp_add_anc.vcf > snp_add_anc_nomiss.vcf #16308
grep -v "##" snp_add_anc.vcf | awk '{print NF; exit}' #216
grep -v "#" snp_add_anc_nomiss.vcf | awk '{for(i=10;i<=NF;i++) count[$i]++} END {for(num in count) print num, count[num]}'
grep -v "#" snp_add_anc.vcf | awk 'BEGIN {A_count=0; D_count=0} {if ($3 == "A") {A_count++} else if ($3 == "D") {D_count++}} END {print "A_count:", A_count, "D_count:", D_count}'  #A_count: 18 D_count: 16559
grep -v "#" snp_add_anc_nomiss.vcf | awk 'BEGIN {A_count=0; D_count=0} {if ($3 == "A") {A_count++} else if ($3 == "D") {D_count++}} END {print "A_count:", A_count, "D_count:", D_count}'  #A_count: 16 D_count: 16232
####derive allele count
bcftools query -l snp_add_anc_nomiss.vcf > xx
    grep -v "#" snp_add_anc_nomiss.vcf | \
    awk 'BEGIN {OFS="\t"} {
        if ($3 == "A") {
            for (i=10; i<=NF; i++) {
                if ($i == "1/1") $i = 2;
                else if ($i == "0/1" || $i == "1/0") $i = 1;
                else if ($i == "0/0") $i = 0;
            }
        } else if ($3 == "D") {
            for (i=10; i<=NF; i++) {
                if ($i == "1/1") $i = 0;
                else if ($i == "0/1" || $i == "1/0") $i = 1;
                else if ($i == "0/0") $i = 2;
            }
        }
        for (i=10; i<=NF; i++) {
            printf "%s%s", $i, (i==NF ? "\n" : OFS)
        }
    }' | awk '{ for(i=1; i<=NF; i++) sum[i] += $i } END { for(i=1; i<=NF; i++) print "sample"i"\t", sum[i] }' | paste xx - | sort -k1,1n > fold0.count.txt
```

### SV
```shell
```bash
##/data/xuebo/pepo/pop_genetics/01_tree/outgroup/SV 使用的是这个
cp SV_nomiss.vcf /data/xuebo/pepo/load/SV
grep -v "#" SV_nomiss.vcf | awk '{print NF; exit}' #有219列,watermelon和Cargyrosperma分别在217和218列
grep -v "#"  SV_nomiss.vcf | awk 'BEGIN {OFS=FS="\t"} {
    if ($217 == "1/1" && $218 == "1/1") {
        $3 = "D"
        print
    } else if ($217 == "0/0" && $218 == "0/0") {
        $3 = "A"
        print
    }
}' > SV_addout #185864
awk 'BEGIN {OFS=FS="\t"} { NF = NF - 3; print}' SV_addout > sv
grep "#" SV_nomiss.vcf > header.txt #手动把三个outgroup去掉
cat header.txt sv > sv_add_anc.vcf 
grep -v "##" sv_add_anc.vcf | awk '{print NF; exit}' #216
grep -v "#" sv_add_anc.vcf | awk '{for(i=10;i<=NF;i++) count[$i]++} END {for(num in count) print num, count[num]}'
grep -v "#" sv_add_anc.vcf | awk 'BEGIN {A_count=0; D_count=0} {if ($3 == "A") {A_count++} else if ($3 == "D") {D_count++}} END {print "A_count:", A_count, "D_count:", D_count}'  #A_count: 121138 D_count: 64726
####derive allele count
bcftools query -l sv_add_anc.vcf > xx
    grep -v "#" sv_add_anc.vcf | \
    awk 'BEGIN {OFS="\t"} {
        if ($3 == "A") {
            for (i=10; i<=NF; i++) {
                if ($i == "1/1") $i = 2;
                else if ($i == "0/1" || $i == "1/0") $i = 1;
                else if ($i == "0/0") $i = 0;
            }
        } else if ($3 == "D") {
            for (i=10; i<=NF; i++) {
                if ($i == "1/1") $i = 0;
                else if ($i == "0/1" || $i == "1/0") $i = 1;
                else if ($i == "0/0") $i = 2;
            }
        }
        for (i=10; i<=NF; i++) {
            printf "%s%s", $i, (i==NF ? "\n" : OFS)
        }
    }' | awk '{ for(i=1; i<=NF; i++) sum[i] += $i } END { for(i=1; i<=NF; i++) print "sample"i"\t", sum[i] }' | paste xx - | sort -k1,1n > sv.count.txt

paste ./SNP_fold4/fold4.count.txt ./SNP_fold0/fold0.count.txt ./SV/sv.count.txt > load207.txt
```

### SNP_adaptive_genotype
```shell
```bash
cp /data/xuebo/pepo/GEA/02_Bayenv/site/top1_SNP/SNP_temp.vcf ./
cp /data/xuebo/pepo/GEA/02_Bayenv/site/top1_SNP/SNP_prec.vcf ./
/data/xuebo/software/vcftools/src/perl/vcf-concat SNP_temp.vcf SNP_prec.vcf > xxxx.vcf
bcftools sort xxxx.vcf -o SNP_adaptive.vcf #6000
bgzip -c SNP_adaptive.vcf > SNP_adaptive.vcf.gz
bcftools annotate -x ^FORMAT/GT -O z -o SNP_adaptive2.vcf.gz SNP_adaptive.vcf.gz
gunzip SNP_adaptive2.vcf.gz
bcftools query -l SNP_adaptive2.vcf > xx
    grep -v "#" SNP_adaptive2.vcf | \
    awk 'BEGIN {OFS="\t"}
     {
         for (i=10; i<=NF; i++) {
             if ($i == "1/1") $i = 2;
             else if ($i == "0/1" || $i == "1/0") $i = 1;
             else if ($i == "0/0") $i = 0;
             else $i = 0;  # 处理./.或缺失值
         }
         for (i=10; i<=NF; i++) {
             sum[i] += $i;
         }
     }
     END {
         for (i=10; i<=NF; i++) {
             print sum[i];
         }
     }' | paste xx - | sort -k1,1n > SNP_adaptive.count.txt
```

### SV_adaptive_genotype
```shell
```bash
cp /data/xuebo/pepo/GEA/02_Bayenv/site/top1_SV/SV_temp.vcf ./
cp /data/xuebo/pepo/GEA/02_Bayenv/site/top1_SV/SV_prec.vcf ./
/data/xuebo/software/vcftools/src/perl/vcf-concat SV_temp.vcf SV_prec.vcf > xxxx.vcf
bcftools sort xxxx.vcf -o SNP_adaptive.vcf #7999
bgzip -c SNP_adaptive.vcf > SNP_adaptive.vcf.gz
bcftools annotate -x ^FORMAT/GT -O z -o SNP_adaptive2.vcf.gz SNP_adaptive.vcf.gz
gunzip SNP_adaptive2.vcf.gz
bcftools query -l SNP_adaptive2.vcf > xx
    grep -v "#" SNP_adaptive2.vcf | \
    awk 'BEGIN {OFS="\t"}
     {
         for (i=10; i<=NF; i++) {
             if ($i == "1/1") $i = 2;
             else if ($i == "0/1" || $i == "1/0") $i = 1;
             else if ($i == "0/0") $i = 0;
             else $i = 0;  # 处理./.或缺失值
         }
         for (i=10; i<=NF; i++) {
             sum[i] += $i;
         }
     }
     END {
         for (i=10; i<=NF; i++) {
             print sum[i];
         }
     }' | paste xx - | sort -k1,1n > SV_adaptive.count.txt

paste ./SNP_adaptive_genotype/SNP_adaptive.count.txt ./SV_adaptive_genotype/SV_adaptive.count.txt > adaptive_genotype207.txt
```



# Suqash from Changlin Wang
文件路径：/data/xuebo/squash
原始数据：/data/xuebo/squash/ori_data
## tree
```shell
```bash
##这个画树，还是使用
pepo.ffds.list #这个是使用C39的4D_site
cp /data/xuebo/pepo/pop_genetics/01_tree/4D_site/pepo.ffds.list ./
cut -f1,2  pepo.ffds.list  | sort -k1,1 -k2,2n > ffds_cut12.list
awk '{print > "4d."$1".txt"}' ffds_cut12.list
sed 's/C39_Chr01/chr1/; s/C39_Chr02/chr2/; s/C39_Chr03/chr3/; s/C39_Chr04/chr4/; s/C39_Chr05/chr5/; s/C39_Chr06/chr6/; s/C39_Chr07/chr7/; s/C39_Chr08/chr8/; s/C39_Chr09/chr9/; s/C39_Chr10/chr10/; s/C39_Chr11/chr11/; s/C39_Chr12/chr12/; s/C39_Chr13/chr13/; s/C39_Chr14/chr14/; s/C39_Chr15/chr15/; s/C39_Chr16/chr16/; s/C39_Chr17/chr17/; s/C39_Chr18/chr18/; s/C39_Chr19/chr19/; s/C39_Chr20/chr20/; s/C39_Chr00/chr0/' squash_SNP_final.vcf > squash_SNP_final2.vcf
nohup bcftools view -R ffds_cut12.list squash_SNP_final2.vcf.gz > fold4_site_SNP_gatk.vcf & ##112894
bgzip -c fold4_site_SNP_gatk.vcf > fold4_site_SNP_gatk.vcf.gz
tabix fold4_site_SNP_gatk.vcf.gz
###加上3个西瓜的外类群
cp /data/xuebo/pepo/pop_genetics/01_tree/outgroup2/SNP_gatk/pepo_add_out3.vcf.gz ./
grep -v "#" fold4_site_SNP_gatk.vcf | cut -f1,2 | sort -k1.4n -k2,2n > ffds_in_vcf.list
nohup bcftools view -R ffds_in_vcf.list pepo_add_out3.vcf.gz > addout3_4D.vcf &
bcftools view -v snps -m2 -M2 -O z -o addout3_4D_2.vcf.gz addout3_4D.vcf
gunzip addout3_4D_2.vcf.gz
bcftools sort addout3_4D_2.vcf -Oz -o addout3_4D_2.vcf.gz
tabix -p vcf addout3_4D_2.vcf.gz
#提取两个 VCF 文件中 REF 和 ALT 都完全相同的变异位点
zcat addout3_4D_2.vcf.gz | awk '!/^#/ {key=$1":"$2; a[key]=$4":"$5} END{for (k in a) print k,a[k]}' > f1.txt
zcat /data/xuebo/pepo/pop_genetics/01_tree/no_outgroup/SNP_gatk/fold4_site_SNP_gatk.vcf.gz | awk '!/^#/ {key=$1":"$2; b[key]=$4":"$5} END{for (k in b) print k,b[k]}' > f2.txt
# 找出 REF/ALT 完全一致的位点
comm -12 <(sort f1.txt) <(sort f2.txt) > shared_variants.txt
awk -F '[: ]' '{print $1"\t"$2}' shared_variants.txt | sort -k1.4n -k2,2n > shared_pos.txt
bcftools view -R shared_pos.txt fold4_site_SNP_gatk.vcf.gz -Oz -o gatk_shared.vcf.gz
bcftools view -R shared_pos.txt addout3_4D_2.vcf.gz -Oz -o addout_shared.vcf.gz
tabix gatk_shared.vcf.gz
tabix addout_shared.vcf.gz
bcftools merge gatk_shared.vcf.gz addout_shared.vcf.gz -Oz -o merged_shared_4D_site_addout.vcf.gz
####做树
/data/xuebo/squash/tree/built_tree
java -jar /data/xuebo/software/WGSc/excuable/WGS.jar --model vcf --type toFasta --file ../merged_shared_4D_site_addout.vcf.gz  --out fold4_site.fasta
java -jar /data/xuebo/software/javaCode/C23_fasta2phy.jar --file1 fold4_site.fasta --out fold4_site.phy
nohup /data/xuebo/miniconda3/envs/orthofinder/bin/raxmlHPC-PTHREADS-SSE3 -f a -m GTRGAMMA -p 12345 -x 12345 -# 100 -s fold4_site.phy -n squash_fold4.raxml -o SRR8751845,SRR8751846,SRR8751851 -T 100 > nohuptree 2>& 1 &
```

## PCA
```shell
```bash                                                                                                                
/data/xuebo/squash/PCA
cp ../tree/fold4_site_SNP_gatk.vcf ./
vcftools --vcf fold4_site_SNP_gatk.vcf --plink --out squash_pca
#PCA
plink --file squash_pca --pca 100 header tabs --chr-set 20 --out squash_pca_PCAout
```

## Admixture
```shell
```bash
/data/xuebo/squash/Admixture
cp ../tree/fold4_site_SNP_gatk.vcf ./
vcftools --vcf fold4_site_SNP_gatk.vcf --plink --out squash_pca
plink --file squash_pca --maf 0.05 --hwe 0.0001 --make-bed --out filter.snp
for i in {1..20};
    do admixture --cv filter.snp.bed $i >> log.txt
done
grep CV log.txt | awk -F ':' '{print NR"\t"$2}' | sed '1i\K\tCV_error' >> CV_for_plot.txt
cut -f1 squash_pca.ped > Admixture_name.txt
for i in {2..7};
do
    paste Admixture_name.txt filter.snp.$i.Q | tr ' ' ',' | tr '\t' ','  > Admixture_4dTv.$i.txt
done
```

## Pangenie
##原始的数据在132:/data/xuebo/squash_Changlin, 之后会删掉的，太大了 
```shell
```bash
cat ../taxtlist.txt | while read line
do
(
java -jar /data/feizj/RNA-Seq/tool/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
    -threads 10 \
    /data/feizj/${line}_R1.fastq.gz \
    /data/feizj/${line}_R2.fastq.gz \
    ${line}_R1_paired.fq.gz \
    ${line}_R1_unpaired.fq.gz \
    ${line}_R2_paired.fq.gz \
    ${line}_R2_unpaired.fq.gz \
    ILLUMINACLIP:/data/feizj/RNA-Seq/tool/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10:2:True \
    LEADING:3 TRAILING:3 MINLEN:36
) &
while [ $(jobs -r | wc -l) -ge 10 ]
do
    sleep 2
done
done
wait
####
for i in *_R*_paired.fq.gz  
do
    num=$(echo $i | cut -d "_" -f1)

    if [ $num -gt 150 ]
    then
        mv $i ../
    fi
done
#########
原始数据在132:/data/xuebo/squash_Changlin/Pangenie
############run PanGenie
132:/data/xuebo/squash_Changlin/Pangenie
#!/bin/bash
MAX_JOBS=4  # 设置每批次并发的最大任务数
N=0         # 当前并发计数器
cat list113.txt | while read line
do
    {
        echo "Processing sample: $line"
        # 解压 reads
        bgzip -@ 50 -cd ../${line}_R1_paired.fq.gz > ${line}_R1.fastq
        bgzip -@ 50 -cd ../${line}_R2_paired.fq.gz > ${line}_R2.fastq
        # 合并 paired reads
        perl /data/xuebo/cucumber/graph/tools/interleave_pairRd.pl ${line}_R1.fastq ${line}_R2.fastq > ${line}_all.fastq
        # 运行 PanGenie
        /data/xuebo/software/pangenie/build/src/PanGenie \
            -f /data/xuebo/squash_Changlin/Pangenie/index/pangenie16.index \
            -i /data/xuebo/squash_Changlin/Pangenie/${line}_all.fastq \
            -s ${line} -j 30 -t 30 \
            -o /data/xuebo/squash_Changlin/Pangenie/output_16haplo_squashwang/${line}
        # 删除中间文件
        rm -f ${line}_R1.fastq ${line}_R2.fastq ${line}_all.fastq
        echo "Finished sample: $line"
    } &  # 放入后台执行
    ((N++))
    if [[ $N -ge $MAX_JOBS ]]; then
        wait    # 等待当前这批任务完成
        N=0     # 重置计数器
    fi
done
wait  # 确保最后一批样本完成
nohup bash getpengenie.sh > nohup_pengenie 2>& 1 &
####
#!/bin/bash
MAX_JOBS=4  # 设置每批次并发的最大任务数
N=0         # 当前并发计数器
cat taxtlist55.txt | while read line
do
    {
        echo "Processing sample: $line"
        # 解压 reads
        bgzip -@ 50 -cd ../clean/${line}_R1_paired.fq.gz > ${line}_R1.fastq
        bgzip -@ 50 -cd ../clean/${line}_R2_paired.fq.gz > ${line}_R2.fastq
        # 合并 paired reads
        perl /data/xuebo/cucumber/graph/tools/interleave_pairRd.pl ${line}_R1.fastq ${line}_R2.fastq > ${line}_all.fastq
        # 运行 PanGenie
        /data/xuebo/software/pangenie/build/src/PanGenie \
            -f /data/xuebo/squash_Changlin/Pangenie/index/pangenie16.index \
            -i /data/xuebo/squash_Changlin/Pangenie/${line}_all.fastq \
            -s ${line} -j 30 -t 30 \
            -o /data/xuebo/squash_Changlin/Pangenie/output_16haplo_squashwang/${line}
        # 删除中间文件
        rm -f ${line}_R1.fastq ${line}_R2.fastq ${line}_all.fastq
        echo "Finished sample: $line"
    } &  # 放入后台执行
    ((N++))
    if [[ $N -ge $MAX_JOBS ]]; then
        wait    # 等待当前这批任务完成
        N=0     # 重置计数器
    fi
done
wait  # 确保最后一批样本完成
nohup bash getpengenie3.sh > nohup_pengenie3 2>& 1 &
####
cat taxtlist.txt | while read line
do
    bgzip -@100 -c ${line}_genotyping.vcf > ${line}_genotyping.vcf.gz
    tabix ${line}_genotyping.vcf.gz
done

/data/xuebo/squash_Changlin/Pangenie/output_16haplo_squashwang
nohup bcftools merge *_genotyping.vcf.gz -o PanGenie_all.vcf & ##这个文件的行数和跑Pangenie的graph的长度一样，5303685
##这个VCF文件太大了，把这个VCF文件简化一下
grep "#" PanGenie_all.vcf > header.txt
grep -v "#" PanGenie_all.vcf | awk 'BEGIN {OFS="\t"} {for(i=1; i<=NF; i++) sub(/:.*/, "", $i)}1' > temp
cat header.txt temp  > PanGenie_all_simplified.vcf ##这个PanGenie得到的所有的variance信息
bcftools norm -m -any  PanGenie_all_simplified.vcf >  PanGenie_all_simplified_bi.vcf ##转成bi-allele
bcftools view PanGenie_all_simplified_bi.vcf | bcftools filter --include 'strlen(REF)<strlen(ALT)' | bcftools view -H > ins & #662061 这个是所有的，后面要挑出来>=20bp的
bcftools view PanGenie_all_simplified_bi.vcf | bcftools filter --include 'strlen(REF)>strlen(ALT)' | bcftools view -H > del & #612689

#########>=20bp SV
#这个文件的问题是缺失表示的是.，现在要把.换成./.
grep -v "#" ins | awk '{if (length($5)-length($4) > 19) print}' | awk 'BEGIN {OFS="\t"} { for (i=10; i<=NF; i++) { if ($i == ".") $i = "./." } print }'> temp
grep "#" PanGenie_all_simplified_bi.vcf > header2.txt
cat header2.txt temp > Insersion_SV20_207_bi.vcf #131373
grep -v "#" del | awk '{if (length($4)-length($5) > 19) print}' | awk 'BEGIN {OFS="\t"} { for (i=10; i<=NF; i++) { if ($i == ".") $i = "./." } print }' > temp
cat header2.txt temp > Deletion_SV20_207_bi.vcf #98058
grep -v "#" snp_mnp | awk '{if (length($4) == 1) print}' | awk 'BEGIN {OFS="\t"} { for (i=10; i<=NF; i++) { if ($i == ".") $i = "./." } print }' > temp
cat header2.txt temp > SNP_207_bi.vcf #3947495
#############################################生成Bi的SV和SNP，带|的，phasing
grep -v "#" Insersion_SV20_207_bi.vcf | awk '{for(i=10;i<=NF;i++) count[$i]++} END {for(num in count) print num, count[num]}'
1/0 129388
1/1 5735606
0/0 20194335
./. 88318
0/1 1046564
##VCF文件排序
(grep "^##" Insersion_SV20_207_bi.vcf
 grep "^#CHROM" Insersion_SV20_207_bi.vcf
 for chr in chr{1..20}; do
     awk -v c=$chr '$1==c' Insersion_SV20_207_bi.vcf
 done) > Insersion_sorted.vcf
(grep "^##" Deletion_SV20_207_bi.vcf
 grep "^#CHROM" Deletion_SV20_207_bi.vcf
 for chr in chr{1..20}; do
     awk -v c=$chr '$1==c' Deletion_SV20_207_bi.vcf
 done) > Deletion_sorted.vcf
 (grep "^##" SNP_207_bi.vcf
 grep "^#CHROM" SNP_207_bi.vcf
 for chr in chr{1..20}; do
     awk -v c=$chr '$1==c' SNP_207_bi.vcf
 done) > SNP_sorted.vcf
##
grep -v "#" Insersion_sorted.vcf | awk 'BEGIN {OFS="\t"} { for (i=10; i<=NF; i++) { if ($i == "0/0") $i = "0|0"; else if ($i == "1/1") $i = "1|1"; else if ($i == "./.") $i = ".|." ;else if ($i == "0/1") $i = "0|0"; else if ($i == "1/0") $i = "0|0"} print }' | cat header2.txt - > Insersion_SV20_207_phasing.vcf
grep -v "#" Deletion_sorted.vcf | awk 'BEGIN {OFS="\t"} { for (i=10; i<=NF; i++) { if ($i == "0/0") $i = "0|0"; else if ($i == "1/1") $i = "1|1"; else if ($i == "./.") $i = ".|." ;else if ($i == "0/1") $i = "0|0"; else if ($i == "1/0") $i = "0|0"} print }' | cat header2.txt - > Deletion_SV20_207_phasing.vcf
grep -v "#" SNP_sorted.vcf | awk 'BEGIN {OFS="\t"} { for (i=10; i<=NF; i++) { if ($i == "0/0") $i = "0|0"; else if ($i == "1/1") $i = "1|1"; else if ($i == "./.") $i = ".|." ;else if ($i == "0/1") $i = "0|0"; else if ($i == "1/0") $i = "0|0"} print }' | cat header2.txt - > SNP_207_phasing.vcf











