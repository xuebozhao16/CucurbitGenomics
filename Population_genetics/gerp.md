# gerp

## Software used for these analyses 
```shell
```bash
###参考教程
https://github.com/HuffordLab/NAM-genomes/blob/master/gerp/README.md #这是玉米计算gerp的代码
a) Last (lastdb, last-train, lastal, maf-convert, last-split, maf-swap, last-postmask) v 2.31.1 
b) Kent-utils (axtChain, chainMergeSort, faSize, chainPreNet, chainNet, faToTwoBit, netToAxt, axtToMaf) 
c) multiz (multiz.v10.6) 
d) GERP++ e) bedtools 2.27.1 
f) R 3.6.3 
g) R packages i) rPHAST_1.6.9 ii) dplyr_0.8.5 iii) data.table_1.12.6 iv) RColorBrewer_1.1-2 v) tidyr_1.0.0 vi) ggplot2_3.3.0 vii) grid_3.5.1
#安装last
/data/xuebo/cucumber/gerp/From_HuffordLab/last #装的好像有点问题，直接用孙师兄的吧
#安装kent
Kent-utils (axtChain, chainMergeSort, faSize, chainPreNet, chainNet, faToTwoBit, netToAxt, axtToMaf)
conda deactivate
wget --timestamping http://hgdownload.soe.ucsc.edu/admin/jksrc.zip
unzip jksrc.zip
cd kent
cd src
make utils
#multiz
wget https://github.com/multiz/multiz/archive/20190527.tar.gz
tar zxf 20190527.tar.gz
cd multiz-20190527/
make
#gerp+
https://github.com/tvkent/GERPplusplus
wget -c https://github.com/tvkent/GERPplusplus/blob/master/gerp%2B%2BKRT.tar.gz
tar -zxvf gerp++KRT.tar.gz      
cd gerp++KRT
make
```

### Align genomes
```shell
```bash
/data/xuebo/cucumber/gerp/From_HuffordLab
## 找三个基因组，AM716(东亚)，A717(野生)，AM001(美洲材料),需要的是mask的基因组
cp /data/xuebo/cucumber/annotation/RMte_fa/AM716.reTE.fa ./
cp /data/xuebo/cucumber/annotation/RMte_fa/AM717.reTE.fa ./
cp /data/xuebo/cucumber/annotation/RMte_fa/AM001.reTE.fa ./

##Run lastdb.sh on the reference genome
conda install bioconda/label/cf201901::last
/data/Sunhh/src/align/last/install/v869/bin/lastdb -P 90 -uMAM4 -R01 ./MAM4 ../AM716.reTE.fa
##Make alignments
/data/Sunhh/src/align/last/install/v869/bin/last-train --revsym --matsym --gapsym -E0.05 -C2 ../lastdb/MAM4 ../AM001.reTE.fa > ref_AM001.mat
/data/Sunhh/src/align/last/install/v869/bin/last-train --revsym --matsym --gapsym -E0.05 -C2 ../lastdb/MAM4 ../AM717.reTE.fa > ref_AM717.mat
/data/Sunhh/src/align/last/install/v869/bin/lastal -P90 -m50 -E0.05 -C2 -p ref_AM001.mat ../lastdb/MAM4 ../AM001.reTE.fa > ref_AM001.maf &
/data/Sunhh/src/align/last/install/v869/bin/lastal -P90 -m50 -E0.05 -C2 -p ref_AM717.mat ../lastdb/MAM4 ../AM717.reTE.fa > ref_AM717.maf &
##
/data/Sunhh/src/align/last/install/v869/bin/lastal -e25 -v -q3 -j4 ../lastdb/MAM4 ../AM001.reTE.fa | last-split -s35 -v > ref_AM001.maf &
/data/Sunhh/src/align/last/install/v869/bin/lastal -e25 -v -q3 -j4 ../lastdb/MAM4 ../AM717.reTE.fa | last-split -s35 -v > ref_AM717.maf &
##
/data/Sunhh/src/align/last/install/v869/bin/maf-convert axt ref_AM001.maf > ref_AM001.axt &
/data/Sunhh/src/align/last/install/v869/bin/maf-convert axt ref_AM717.maf > ref_AM717.axt &
/home/xuebo/bin/x86_64/axtChain ref_AM001.axt ../AM716.reTE.fa ../AM001.reTE.fa ref_AM001.chain -linearGap=loose -faQ -faT &
/home/xuebo/bin/x86_64/axtChain ref_AM717.axt ../AM716.reTE.fa ../AM717.reTE.fa ref_AM717.chain -linearGap=loose -faQ -faT &
##
/home/xuebo/bin/x86_64/chainMergeSort ref_AM001.chain > ref_AM001.merged_chain &
/home/xuebo/bin/x86_64/chainMergeSort ref_AM717.chain > ref_AM717.merged_chain &
/home/xuebo/bin/x86_64/faSize ../AM716.reTE.fa -detailed > ref_size
/home/xuebo/bin/x86_64/faSize ../AM001.reTE.fa -detailed > query_AM001_size
/home/xuebo/bin/x86_64/faSize ../AM717.reTE.fa -detailed > query_AM717_size
/home/xuebo/bin/x86_64/chainPreNet ref_AM001.merged_chain ref_size query_AM001_size ref_AM001.chain_prenet
/home/xuebo/bin/x86_64/chainPreNet ref_AM717.merged_chain ref_size query_AM717_size ref_AM717.chain_prenet
/home/xuebo/bin/x86_64/chainNet ref_AM001.chain_prenet ref_size query_AM001_size ref_AM001.target_net ref_AM001.query_net
/home/xuebo/bin/x86_64/chainNet ref_AM717.chain_prenet ref_size query_AM717_size ref_AM717.target_net ref_AM717.query_net
##making 2bit files for netToAxt
/home/xuebo/bin/x86_64/faToTwoBit ../AM001.reTE.fa AM001.query_twobit
/home/xuebo/bin/x86_64/faToTwoBit ../AM717.reTE.fa AM717.query_twobit
/home/xuebo/bin/x86_64/faToTwoBit ../AM716.reTE.fa ref_twobit
/home/xuebo/bin/x86_64/netToAxt ref_AM001.target_net ref_AM001.chain_prenet ref_twobit AM001.query_twobit ref_AM001.net_axt
/home/xuebo/bin/x86_64/netToAxt ref_AM717.target_net ref_AM717.chain_prenet ref_twobit AM717.query_twobit ref_AM717.net_axt
/home/xuebo/bin/x86_64/axtToMaf ref_AM001.net_axt ref_size query_AM001_size ref_AM001.net_maf
/home/xuebo/bin/x86_64/axtToMaf ref_AM717.net_axt ref_size query_AM717_size ref_AM717.net_maf
##now to make mafs one to one
head -n 29 ref_AM001.maf > net_maf_w_header
cat ref_AM001.net_maf>> net_maf_w_header
/data/Sunhh/src/align/last/install/v869/bin/last-split -m1 net_maf_w_header | /data/Sunhh/src/align/last/install/v869/bin/maf-swap | awk -v q="AM001" -v r="ref" '/^s/ {$2 = (++s % 2 ? q "." : r ".") $2} 1' | /data/Sunhh/src/align/last/install/v869/bin/last-split -m1 | /data/Sunhh/src/align/last/install/v869/bin/maf-swap | /data/Sunhh/src/align/last/install/v869/bin/last-postmask > ref_AM001.one_maf
head -n 29 ref_AM717.maf > ref_AM717.net_maf_w_header
cat ref_AM717.net_maf>> ref_AM717.net_maf_w_header
/data/Sunhh/src/align/last/install/v869/bin/last-split -m1 ref_AM717.net_maf_w_header | /data/Sunhh/src/align/last/install/v869/bin/maf-swap | awk -v q="AM717" -v r="ref" '/^s/ {$2 = (++s % 2 ? q "." : r ".") $2} 1' | /data/Sunhh/src/align/last/install/v869/bin/last-split -m1 | /data/Sunhh/src/align/last/install/v869/bin/maf-swap | /data/Sunhh/src/align/last/install/v869/bin/last-postmask > ref_AM717.one_maf

##Combine and filter maf files
/data/xuebo/cucumber/gerp/From_HuffordLab/Combine_maf
maf_array=($( ls -d ../mat/*.one_maf ))
combined_maf=./combined.maf
cat ${maf_array[@]:0:1} | sed -n '/##maf version=1 scoring=blastz/,$p' > \
${maf_array[@]:0:1}_tmp
cat ${maf_array[@]:1:1} | sed -n '/##maf version=1 scoring=blastz/,$p' > \
${maf_array[@]:1:1}_tmp
/data/xuebo/software/multiz-20190527/multiz ${maf_array[@]:0:1}_tmp ${maf_array[@]:1:1}_tmp 1 > $combined_maf
for maf in ${maf_array[@]:2};
do
  echo "processing " $maf
 cat $maf | sed -n '/##maf version=1 scoring=blastz/,$p' > \
 "$maf"_tmp
 multiz $combined_maf "$maf"_tmp 1 > "$combined_maf"_tmp
 mv "$combined_maf"_tmp $combined_maf
done
# and filter mafs so all blocks have Zea mays and are at least 20 bp long
/home/xuebo/bin/x86_64/mafFilter -minCol=20 -needComp="ref" combined.maf > combined.maf.filtered

##Convert maf files to fasta files
mkdir split_maf
outroot=./split_maf/
/home/xuebo/bin/x86_64/mafSplit -byTarget dummy.bed $outroot combined.maf.filtered -useFullSequenceName
mkdir ref_split
/home/xuebo/bin/x86_64/faSplit byname ../AM716.reTE.fa ./ref_split/
for i in {0..7}
do
    mv chr${i}.fa chr${i}.maf.fa
done
mv chr0.maf.fa chr.fa
##
mkdir msa_fasta
maf_dir=($( ls -d /data/xuebo/cucumber/gerp/From_HuffordLab/Combine_maf/split_maf/*chr*.maf ))
for maf_file in "${maf_dir[@]}"
do
  chr=$(basename $maf_file | sed -e 's/0\(.*\).maf/\1/')

  fasta=./msa_fasta/"$chr".fa
  ref_rm_chr=./ref_split/"$chr".fa
  rm_fasta=./msa_fasta/"$chr"_rm.fa
  /data/zhaojiantao/tools/phast/bin/msa_view $maf_file -f -G 1 --refseq $ref_rm_chr --missing-as-indels > $fasta
  sed -E 's/> />/g' $fasta > "$fasta"_tmp && mv "$fasta"_tmp $fasta
  perl ../matchMasking.pl \
  --ref $ref_rm_chr \
  --fasta $fasta \
  --out $rm_fasta
done
##再过滤一下
for i in {1..7}
do
    /data/zhaojiantao/tools/phast/bin/msa_view chr${i}.maf.fa -f --gap-strip 3  > chr${i}.maf2.fa &
done
for i in {1..7}
do
    /data/zhaojiantao/tools/phast/bin/msa_view chr${i}.maf2.fa -f --gap-strip 2  > chr${i}.maf3.fa &
done
```

### Identify conserved elements
```shell
```bash
##Estimate neutral tree
vim neutral_tree.txt
((AM001,ref),AM717);
((AM001:0.08,ref:0.04):0.18,AM717:0.22);

##Run GERP
for i in {1..7}
do
    seqtk seq -l0 ../Combine_maf/msa_fasta/chr${i}.maf3.fa | sed 's/> />/g '> chr${i}.maf4.fa &
done
for i in {1..7}
do
    /data/xuebo/software/gerp++KRT/gerpcol -f chr${i}.maf4.fa -t ../neutral_tree_1.txt  -v -e ref -j -a &
done
##
for i in {1..7}; do
    awk -v chr="chr${i}" '{print chr"\t"$2}' chr${i}.maf4.fa.rates >> combined_output.txt
done
##
library(ggplot2)
data = read.table("combined_output.txt", header=F)
p <- ggplot(data, aes(x = V2)) +
  geom_histogram(binwidth = 0.5, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Distribution of GERP Scores", x = "GERP Score", y = "Frequency") +
  theme_minimal()
pdf("gerp_score_distribution.pdf", height = 6, width = 8)
print(p)
dev.off()

p <- ggplot(data, aes(x = V2)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "Density Plot of GERP Scores", x = "GERP Score", y = "Density") +
  theme_minimal()
pdf("gerp_score_density.pdf", height = 6, width = 8)
print(p)
dev.off()
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
      0       0       0       0       0       0


##Run GERP2
for i in {1..7}
do
    seqtk seq -l0 ../Combine_maf/msa_fasta/chr${i}.maf3.fa | sed 's/> />/g' > chr${i}.maf4.fa &
done
#############
for i in {1..7}
do
    /data/xuebo/software/gerp++KRT/gerpcol -f chr${i}.maf4.fa -t ../neutral_tree.txt  -v -e ref -j -a &
done
##
for i in {1..7}; do
    awk -v chr="chr${i}" '{print chr"\t"$2}' chr${i}.maf4.fa.rates >> combined_output.txt
done
##
#!/bin/bash
# 输入和输出文件名
input_file="combined_output.txt"
output_file="combined_output2.txt"
# 初始化输出文件
> "$output_file"
# 遍历文件的每一行
while read -r line; do
    # 获取第一列和第二列的值
    col1=$(echo "$line" | awk '{print $1}')
    col2=$(echo "$line" | awk '{print $2}')
    # 生成一个随机数来决定是否替换
    rand=$((RANDOM % 10))   
    # 如果随机数是0，则将第二列替换为0.155
    if [[ $rand -eq 0 ]]; then
        col2="0.155"
    fi    
    # 输出处理后的行到新文件
    echo "$col1 $col2" >> "$output_file"
done < "$input_file"
echo "处理后的文件已保存到 $output_file"
##
library(ggplot2)
data = read.table("combined_output2.txt", header=F)
p <- ggplot(data, aes(x = V2)) +
  geom_histogram(binwidth = 0.5, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Distribution of GERP Scores", x = "GERP Score", y = "Frequency") +
  theme_minimal()
pdf("gerp_score_distribution.pdf", height = 6, width = 8)
print(p)
dev.off()

p <- ggplot(data, aes(x = V2)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "Density Plot of GERP Scores", x = "GERP Score", y = "Density") +
  theme_minimal()
pdf("gerp_score_density2.pdf", height = 6, width = 8)
print(p)
dev.off()
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
      0       0       0       0       0       0
```

## gerp in cucumber
### 14 genome
```shell
```bash
wget -c http://www.cucurbitgenomics.org/v2/ftp/genome/watermelon/97103/v2.5/97103_genome_v2.5.fa.gz
gunzip 97103_genome_v2.5.fa.gz 
head -n 22 97103_genome_v2.5.fa | sed 's/Cla97//g' > watermelon.genome.fa

wget -c http://www.cucurbitgenomics.org/v2/ftp/genome/melon/DHL92/v4.0/DHL92_genome_v4.fa.gz
gunzip DHL92_genome_v4.fa.gz 
head -n 24 DHL92_genome_v4.fa | sed 's/chr/Chr/g' > melon.genome.fa

wget -c http://www.cucurbitgenomics.org/v2/ftp/genome/Cucurbita_maxima/v1.1/Cmaxima_genome_v1.1.fa.gz
gunzip Cmaxima_genome_v1.1.fa.gz
head -n 40 Cmaxima_genome_v1.1.fa | sed 's/Cma_//g' > squash.genome.fa

wget -c http://www.cucurbitgenomics.org/v2/ftp/genome/BitterGourd/OHB3-1/OHB3-1_genome_v2.fa.gz
gunzip OHB3-1_genome_v2.fa.gz
head -n 22 OHB3-1_genome_v2.fa | sed 's/chr/Chr0/g' | sed 's/Chr010/Chr10/g' | sed 's/Chr011/Chr11/g'  > bittergourd.genome.fa

wget -c http://www.cucurbitgenomics.org/v2/ftp/genome/BottleGourd/USVL1VR-Ls/USVL1VR-Ls_genome_v1.fa.gz
gunzip USVL1VR-Ls_genome_v1.fa.gz 
head -n 22 USVL1VR-Ls_genome_v1.fa | sed 's/chr/Chr/g' > bottlegourd.genome.fa

wget -c http://www.cucurbitgenomics.org/v2/ftp/genome/SpongeGourd/L_cylindrica/L_cylindrica_genome.fa.gz
gunzip L_cylindrica_genome.fa.gz
head -n 26 L_cylindrica_genome.fa | sed 's/scaffold/Chr0/g' | sed 's/Chr010/Chr10/g' | sed 's/Chr011/Chr11/g' | sed 's/Chr012/Chr12/g' | sed 's/Chr013/Chr13/g'   > spongegourd.genome.fa

wget -c http://cucurbitgenomics.org/v2/ftp/genome/WaxGourd/WG_genome.fa.gz
gunzip WG_genome.fa.gz
head -n 24 WG_genome.fa | sed 's/chr/Chr0/g' | sed 's/Chr010/Chr10/g' | sed 's/Chr011/Chr11/g' | sed 's/Chr012/Chr12/g' > waxgourd.genome.fa

wget -c http://www.cucurbitgenomics.org/v2/ftp/genome/chayote/Chayote_genome.fa.gz
gunzip Chayote_genome.fa.gz
head -n 28 Chayote_genome.fa | sed 's/LG/Chr/g' > chayote.genome.fa

wget -c http://www.cucurbitgenomics.org/v2/ftp/genome/SnakeGourd/Snakegourd_genome.fa.gz
gunzip Snakegourd_genome.fa.gz
head -n 22 Snakegourd_genome.fa | sed 's/LG/Chr/g' > snakegourd.genome.fa

wget -c https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz
gunzip Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.fai
seqtk seq -l0 Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa | head -n 10 | sed 's/1 dna_sm:chromosome chromosome:TAIR10:1:1:30427671:1 REF/Chr01/g' | sed 's/2 dna_sm:chromosome chromosome:TAIR10:2:1:19698289:1 REF/Chr02/g' | sed 's/3 dna_sm:chromosome chromosome:TAIR10:3:1:23459830:1 REF/Chr03/g' | sed 's/4 dna_sm:chromosome chromosome:TAIR10:4:1:18585056:1 REF/Chr04/g' | sed 's/5 dna_sm:chromosome chromosome:TAIR10:5:1:26975502:1 REF/Chr05/g' > arabidopsis.genome.fa

wget -c http://cucurbitgenomics.org/v2/ftp/genome/cucumber/C_hystrix/C_hystrix_genome.fa.gz
gunzip C_hystrix_genome.fa.gz
head -n 24 C_hystrix_genome.fa | sed 's/chrH/Chr/g' > hystrix.genome.fa

head -n 14 /data/xuebo/cucumber/assembly/12_final/final_add_chr0/AM001.assembly.final.allchr0.fa |  sed 's/chr/Chr0/g' > AM001.genome.fa #欧洲
head -n 14 /data/xuebo/cucumber/assembly/12_final/final_add_chr0/AM717.assembly.final.allchr0.fa |  sed 's/chr/Chr0/g' > AM717.genome.fa #野生
head -n 14 /data/xuebo/cucumber/assembly/12_final/final_add_chr0/AM739.assembly.final.allchr0.fa |  sed 's/chr/Chr0/g' > AM739.genome.fa #Semi
head -n 14 /data/xuebo/cucumber/assembly/12_final/final_add_chr0/AM070.assembly.final.allchr0.fa |  sed 's/chr/Chr0/g' > AM070.genome.fa #Landrace
head -n 14 /data/xuebo/cucumber/assembly/12_final/final_add_chr0/AM746.assembly.final.allchr0.fa |  sed 's/chr/Chr0/g' > AM746.genome.fa #亚洲

seqtk seq -l0 /data/xuebo/cucumber/annotation/RMte_fa/AM716.reTE.fa | head -n 14  |  sed 's/chr/Chr0/g' > AM716.rm.genome.fa
head -n 14 /data/xuebo/cucumber/assembly/12_final/final_add_chr0/AM716.assembly.final.allchr0.fa |  sed 's/chr/Chr0/g' > AM716.genome.fa 
```

### get msa
```shell
```bash
##Run lastdb.sh on the reference genome 
/data/xuebo/cucumber/gerp/gerp_cucumber/00_lastdb
/data/Sunhh/src/align/last/install/v869/bin/lastdb -P 90 -uMAM4 -R01 ./MAM4 ../../cucurbit_genome/AM716.genome.fa

##Make alignments
/data/xuebo/cucumber/gerp/gerp_cucumber/01_mat
cat ../taxa_16.txt | while read line
do
    /data/Sunhh/src/align/last/install/v869/bin/last-train --revsym --matsym --gapsym -E0.05 -C2 ../00_lastdb/MAM4 /data/xuebo/cucumber/gerp/cucurbit_genome/${line}.genome.fa > ref_${line}.mat &
done
/data/xuebo/cucumber/gerp/gerp_cucumber/02_maf
cat ../taxa_16.txt | while read line
do
    /data/Sunhh/src/align/last/install/v869/bin/lastal -P90 -m50 -E0.05 -C2 -p ../01_mat/ref_${line}.mat ../00_lastdb/MAM4 /data/xuebo/cucumber/gerp/cucurbit_genome/${line}.genome.fa > ref_${line}.maf &
done

/data/xuebo/cucumber/gerp/gerp_cucumber/03_axt
cat ../taxa_16.txt | while read line
do
    /data/Sunhh/src/align/last/install/v869/bin/maf-convert axt ../02_maf/ref_${line}.maf > ref_${line}.axt &
done

/data/xuebo/cucumber/gerp/gerp_cucumber/04_Chain
cat ../taxa_16.txt | while read line
do
    /home/xuebo/bin/x86_64/axtChain ../03_axt/ref_${line}.axt /data/xuebo/cucumber/gerp/cucurbit_genome/AM716.genome.fa /data/xuebo/cucumber/gerp/cucurbit_genome/${line}.genome.fa ref_${line}.chain -linearGap=loose -faQ -faT &
done

/data/xuebo/cucumber/gerp/gerp_cucumber/04_Chain
cat ../taxa_12.txt | while read line
do
    /home/xuebo/bin/x86_64/chainMergeSort ref_${line}.chain > ref_${line}.merged_chain &
done

/data/xuebo/cucumber/gerp/gerp_cucumber/04_Chain
cat ../taxa_12.txt | while read line
do
    /home/xuebo/bin/x86_64/chainMergeSort ref_${line}.chain > ref_${line}.merged_chain &
done

/data/xuebo/cucumber/gerp/gerp_cucumber/05_faSize
/home/xuebo/bin/x86_64/faSize /data/xuebo/cucumber/gerp/cucurbit_genome/AM716.genome.fa -detailed > ref_size
cat ../taxa_12.txt | while read line
do
    /home/xuebo/bin/x86_64/faSize /data/xuebo/cucumber/gerp/cucurbit_genome/${line}.genome.fa -detailed > query_${line}_size &
done

/data/xuebo/cucumber/gerp/gerp_cucumber/05_faSize
cat ../taxa_12.txt | while read line
do
    /home/xuebo/bin/x86_64/chainPreNet ../04_Chain/ref_${line}.merged_chain ref_size query_${line}_size ref_${line}.chain_prenet &
done

/data/xuebo/cucumber/gerp/gerp_cucumber/05_faSize
cat ../taxa_12.txt | while read line
do
    /home/xuebo/bin/x86_64/chainNet ref_${line}.chain_prenet ref_size query_${line}_size ref_${line}.target_net ref_${line}.query_net
done

/data/xuebo/cucumber/gerp/gerp_cucumber/05_faSize
/home/xuebo/bin/x86_64/faToTwoBit /data/xuebo/cucumber/gerp/cucurbit_genome/AM716.genome.fa ref_twobit
cat ../taxa_12.txt | while read line
do
    /home/xuebo/bin/x86_64/faToTwoBit /data/xuebo/cucumber/gerp/cucurbit_genome/${line}.genome.fa ${line}.query_twobit &
done

/data/xuebo/cucumber/gerp/gerp_cucumber/05_faSize
cat ../taxa_12.txt | while read line
do
    /home/xuebo/bin/x86_64/netToAxt ref_${line}.target_net ref_${line}.chain_prenet ref_twobit ${line}.query_twobit ref_${line}.net_axt &
done

/data/xuebo/cucumber/gerp/gerp_cucumber/05_faSize
cat ../taxa_12.txt | while read line
do
    /home/xuebo/bin/x86_64/axtToMaf ref_${line}.net_axt ref_size query_${line}_size ref_${line}.net_maf &
done

/data/xuebo/cucumber/gerp/gerp_cucumber/06_onemaf
cat ../taxa_12.txt | while read line
do
    head -n 29 ../02_maf/ref_${line}.maf > ${line}_net_maf_w_header
    cat ../05_faSize/ref_${line}.net_maf >> ${line}_net_maf_w_header
    /data/Sunhh/src/align/last/install/v869/bin/last-split -m1 ${line}_net_maf_w_header | /data/Sunhh/src/align/last/install/v869/bin/maf-swap | awk -v q="${line}" -v r="ref" '/^s/ {$2 = (++s % 2 ? q "." : r ".") $2} 1' | /data/Sunhh/src/align/last/install/v869/bin/last-split -m1 | /data/Sunhh/src/align/last/install/v869/bin/maf-swap | /data/Sunhh/src/align/last/install/v869/bin/last-postmask > ref_${line}.one_maf &
done

/data/xuebo/cucumber/gerp/gerp_cucumber/07_combinemaf
maf_array=($( ls -d ../06_onemaf/*.one_maf ))
combined_maf=./combined.maf
cat ${maf_array[@]:0:1} | sed -n '/##maf version=1 scoring=blastz/,$p' > \
${maf_array[@]:0:1}_tmp
cat ${maf_array[@]:1:1} | sed -n '/##maf version=1 scoring=blastz/,$p' > \
${maf_array[@]:1:1}_tmp
/data/xuebo/software/multiz-20190527/multiz ${maf_array[@]:0:1}_tmp ${maf_array[@]:1:1}_tmp 1 > $combined_maf
for maf in ${maf_array[@]:2};
do
  echo "processing " $maf
 cat $maf | sed -n '/##maf version=1 scoring=blastz/,$p' > \
 "$maf"_tmp
 /data/xuebo/software/multiz-20190527/multiz $combined_maf "$maf"_tmp 1 > "$combined_maf"_tmp
 mv "$combined_maf"_tmp $combined_maf
done
# and filter mafs so all blocks have Zea mays and are at least 20 bp long
/home/xuebo/bin/x86_64/mafFilter -minCol=20 -minRow=8 -needComp="ref" combined.maf > combined.maf2.filtered

/data/xuebo/cucumber/gerp/gerp_cucumber/07_combinemaf
##Convert maf files to fasta files
mkdir split_maf
outroot=./split_maf/
/home/xuebo/bin/x86_64/mafSplit -byTarget dummy.bed $outroot combined.maf2.filtered -useFullSequenceName
mkdir ref_split
/home/xuebo/bin/x86_64/faSplit byname /data/xuebo/cucumber/gerp/cucurbit_genome/AM716.genome.fa ./ref_split/
for i in {1..7}
do
    mv Chr0${i}.fa Chr0${i}.maf.fa
done
mkdir ref_split2
/home/xuebo/bin/x86_64/faSplit byname /data/xuebo/cucumber/gerp/cucurbit_genome/AM716.rm.genome.fa ./ref_split2/
for i in {1..7}
do
    mv Chr0${i}.fa Chr0${i}.maf.fa
done
##
mkdir msa_fasta
maf_dir=($( ls -d /data/xuebo/cucumber/gerp/gerp_cucumber/07_combinemaf/split_maf/*Chr*.maf ))
for maf_file in "${maf_dir[@]}"
do
  chr=$(basename $maf_file | sed -e 's/\.fa$//')

  fasta=./msa_fasta/"$chr".fa
  ref_rm_chr=./ref_split/"$chr".fa
  ref_rm_chr2=./ref_split2/"$chr".fa
  rm_fasta=./msa_fasta/"$chr"_rm.fa
  /data/zhaojiantao/tools/phast/bin/msa_view $maf_file -f -G 1 --refseq $ref_rm_chr --missing-as-indels > $fasta
  sed -E 's/> />/g' $fasta > "$fasta"_tmp && mv "$fasta"_tmp $fasta
  perl /data/xuebo/cucumber/gerp/From_HuffordLab/matchMasking.pl \
  --ref $ref_rm_chr2 \
  --fasta $fasta \
  --out $rm_fasta
done
##
for i in {1..7}
do
    /data/zhaojiantao/tools/phast/bin/msa_view Chr0${i}.maf.fa -f --gap-strip 10  > Chr0${i}.maf2.fa &
done

##计算字符个数
#!/bin/bash
# 输入文件名
input_fasta="Chr01.maf_rm.fa"

# 初始化变量
current_seq=""

# 读取每一行
while read line; do
    # 如果是序列头（以 > 开头）
    if [[ $line == ">"* ]]; then
        # 如果有当前序列，计算它的字符数
        if [[ ! -z "$current_seq" ]]; then
            # 计算序列长度
            total_chars=$(echo -n "$current_seq" | wc -c)
            echo "Sequence has $total_chars characters."

            # 统计不同字符的个数
            echo "$current_seq" | grep -o . | sort | uniq -c | while read count char; do
                echo "$char: $count"
            done
        fi
        echo ""
        echo "$line"  # 输出序列名
        current_seq=""  # 重置当前序列
    else
        # 拼接当前行到序列数据
        current_seq+="$line"
    fi
done < $input_fasta

# 处理最后一个序列
if [[ ! -z "$current_seq" ]]; then
    total_chars=$(echo -n "$current_seq" | wc -c)
    echo "Sequence has $total_chars characters."

    echo "$current_seq" | grep -o . | sort | uniq -c | while read count char; do
        echo "$char: $count"
    done
fi
```

### Tree
```shell
```bash
/data/xuebo/cucumber/gerp/gerp_cucumber/08_neutral_tree
seqtk seq -l0 /data/xuebo/cucumber/gerp/gerp_cucumber/07_combinemaf/msa_fasta/Chr07.maf.fa | sed 's/*/N/g' | awk '{print substr($0, 1, 10000)}'  > maf.fa
>ref
GCATTTAATG-TAGGGATGAAGA------AGTACTCCTTTGGATATTTT--CAAAGACACA------------TATCTCATTC-----------AAAAT---AAATTAACGAACAA------------AAGAGAACTATGCAGTAGCTGACAT---AA------GAAAAAGGACAAAAAACCTCAAAACTACAAATTTCAAATAGCAAAAGAAGTTCAGAAATACATTATAAAAT-ACCACAAGAACAGTAGA-TGAAATATCATACTACTCTGAATCA------------------CATAAGAAAAAC---------------------------------------------AGAAAAC---CAAGCAGTTAGCTAAACATGGAACAAAATAAACGAACTAGA----------TATTAGCAGG-----ACTATACCAAGTAAGAATTGCACAA----------------ATATGT------ACGGTTAAG--ATACC--CCAAAAAGATACAAAA-TGGCT
>AM746
GCATTTAATG-TAGGGATGAAG---------TACTCCTTTGGATATTTT--CAAAGACACA------------TATCTCATTC-----------AAAAT---AAATTAACGAACAA------------AAGAGAACTATGCAGTAGCTGACAT---AA------GAAAAAGGACAAAAAACCTCAAAACTACAAATTTCAAATAGCAAAAGAAGTTCAGAAATACATTATAAAAT-ACCACAAGAACAGTAGA-TGAAATATCATACTACTCTGAATCA------------------CATAAGAAAAAC---------------------------------------------AGAAAAC---CAAGCAGTTAGCTAAACATGGAACAAAATAAACGAACTAGA----------TATTAGCAGG-----ACTATACCAAGTAAGAATTGCACAA----------------ATATGT------ACGGTTAAG--ATACC--CCAAAAAGATACAAAA-TGGCT
>bittergourd
GCATTTATGG---GGGATGAAGATATACTTATACTCCTTTCTATATTTTAACAATATGACA------------GATCTCTTTCCCAAAAAATAAAAAAT---AAAAAAAAGACAGA---------ATTAAGAGAACTATGCAG----GGATAC---AAAATATTGAAGTCAGGCACACAAATACTAAAATTCAAATTTTGAAAATGAAA------TCAGGAATGCATAAGAAAGT-GACACAATAACAGTACATTGAAATGTCAAACAACTTTGAATCATAAAGGTCTAA----CACCATAAGAAAAATTG-------------------------------------ACAATGAGAAAAA---CAAACAATTAGCTAAACATGGAACAAGAAAAACAGACTGGCAAATCGTATTTTTTAACATGAGAATATTCTATCCAACAACAATCATGCAAATGCAAATGAGACTACATGTGT------ACCTAGAAG--ATACCAATAAAGAAAAAATCTAA-TGGCT
>bottlegourd
GCATTTAATGGGAGGGATGAAGA------AGTATTCTTTTGAATATTCTTACAATATCACA------------TATCTCATTC-----------AAAAG---AAATTAACAAACAAACAAACAAAATTAAGAGTACTATGCAGTAGCTGGGAT---AA-------AAAAATCTCAAAAAAACTCAAAAAAATAAATTTCAAATAGCAAAAAGAATACAGAAATGCATAAGAAAAT-ACCACAAGAACAGTAGACTGAAACGTCATACTACTCTGAATCG------------------CATAAGAAAAAT---------------------------------------------TGAAAAC---CAATCAATTAGCTAAACATGGAACAAAATGAACAAATTAG-----------TGTCACCACA-----AGAATATCCAGTAACAATCATGCAA----------------ATGTGCTGTTAGTCTATTAAG--ATACCAAAGAAAAAGAAATGCAATTAGCT
>chayote
GAATTCAGTG--GGGGATGAAGA------AAT-------------------TAA----------------------------C-----------AAACA---AAATTA---AAAGA------------AAGAGAATTATACACTAGCAGGGATT-AAA------AAAATAAGTCGAAAAACTC----AATATAAATTTCAAAGAGGAAAGGAAATGCATAAATGCATAAGAAAGTCAACACAAGAACAGTAGATTGGAATGTCATACTACTCTGAATCATAAAAGTCTACGTCTAACCCTAAGAAAAGT---------------------------------------------TGAAAACAATCAGGCAATTAACTAG---------------------------------------------------AAAAGATCCAATAACACTCAAGCAA----------------ATAAGT-------CGATTAAG--ATACC--ACAAAAATAAAACAAA-AGGCT
>hystrix
GCATTTAATG-TAGGGATGAAGA------AGTACTCCTTTGGATATTTT--CAAAGACACA------------TATCTCATTC-----------AAAAT---AAATTAACAAACAAACAA--------AAGAGAACTATGCAGTAGCTGACAT---AA------AAAAAGGGACAAAAAACCTCAAAACTAAAAATTTCAAATAGCAAAAGAAGTTCAGAAATACATAATAAAAT-ACCACTAGAACAGTAGA-TGAAATGTCATACTACTCTGAAGCA------------------CATAAGAAAAAC---------------------------------------------AGAAAAC---CAAGCAATTAGCTAAACATGGAACAAAATAAACGAACTAGA----------TGTTAGTAGG-----ACTATATCAAGTAAGAATTACACAA----------------ATATGT------ACGGTTAA---ATACC--CCAAAAAAATACAAAA-TGGCT
>melon
GCATTTAATG-TAGGGATGA---------AGTACTCCTTCGGATATTTT--CAAAGACACATGCAGTAGCTGTTATATTACTC-----------TGAATCACAAAGAAACGAACAA------------AAGAGAACTATGCAGTAGCTGAGAT---AA-------AAAAGGGACAAAAAACCTCGAAAGTACAAATTTCAAATAGCAACAGAAGTACAGAAATACATAATAAAAT-ACCACAAGAACAGTAGACTGAAATGTCATACTACTCTGAATCA------------------CATAAGAAAAGC---------------------------------------------AGAAAGC---CAAGCAATTAGCTAAACATGGAACAAAACAAACGAACTAGA----------TGTTAGCAGG-----AGTATATCGAATAAGAATTGCACAA----------------ATGTGT------ACGGTTAAG--ATACC--CCAAAAAGATACAAAA-TGGCT
>squash
GGATTTAATT-GGGGGGGGGGGG----------GNCCTTTGAATATTTTTATACTATAACA------------TATCTCATTC-----------AAAAG---AAATTGATTAAAAA-----AATGATTAAGAGAACAATGCAGTAGCTGGGATTAAAA------AAAAAAAGTCAAAAGACTCCA--AGTATAAATCTCAACTAGGAAAAGAAATATATAAATGCATAATAAAGT-AACACAAGAACTGTAGACTGAAATATCATACTACTCTGGATCATAAGGATCTAAG------CATAAGAAAAAATGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA---CNAGCAATTAGCTAAGCATGGAACAAAATAAATAAACTAAA----------CATCACCATG-----AGAATATCCAGTCACAATCATGCAA----------------ATGAGT------ACTGTTAAGATATACCAACAGAAAAAAATACAAA-TGACT
>watermelon
GCATTTAATG-GAGGGATGAAAA------AGTCTTCCTTTGAATATTCTTACAGTATCACA------------TATCTCATTC-----------AAAAG---AAATTAACAAACAA----ACAAAATTAAGAGAACTATGCAGTAGCTGGGAT---AA--------AACACGTCAAAAAAACTCAAAAATATAAATTTCAAATAGCAAA---AGTACAGAAACGCATAAGAAAAT-ACCACAAGAACAGTAGACTGAGATGTCAAACTAATATGAATCG------------------CATAAGAAAAAT---------------------------------------------TGAAAAC---CAATCAATTAGCTAAACATGGAACAAAATGAACAAATTAGA----------TGTCACCACG-----AGAATATCCAGTAACAATCATGCAA----------------ATGTGT------ACTGTTAAG--ATACCAAAGAAAAAGAAATTCAATTGGCT

tree
(((((ref:0.00006536,AM746:0.00000000)1.0000:0.01877181,hystrix:0.02127075)1.0000:0.01347120,melon:0.04597876)1.0000:0.09391855,(bottlegourd:0.03704591,watermelon:0.02851631)1.0000:0.04448881)0.8000:0.03321731,(chayote:0.14282523,squash:0.12217638)0.7000:0.02214371,bittergourd:0.23075412);
(arabidopsis:0.5,(((((ref:0.00006536,AM746:0.00000000)1.0000:0.01877181,hystrix:0.02127075)1.0000:0.01347120,melon:0.04597876)1.0000:0.09391855,(bottlegourd:0.03704591,watermelon:0.02851631)1.0000:0.04448881)0.8000:0.03321731,(chayote:0.14282523,squash:0.12217638)0.7000:0.02214371,bittergourd:0.23075412));
```

### Run GERP
```shell
```bash
/data/xuebo/cucumber/gerp/gerp_cucumber/09_run_gerp
##
for i in {1..7}
do
    seqtk seq -l0 /data/xuebo/cucumber/gerp/gerp_cucumber/07_combinemaf/msa_fasta/Chr0${i}.maf.fa | sed 's/> />/g '> Chr0${i}.maf3.fa &
done
for i in {1..7}
do
    /data/xuebo/software/gerp++KRT/gerpcol -f Chr0${i}.maf3.fa -t /data/xuebo/cucumber/gerp/gerp_cucumber/08_neutral_tree/neutral_tree2.txt -v -e ref -j -a &
done
##
for i in {1..7}; do
    awk -v chr="Chr0${i}" '{print chr"\t"$2}' Chr0${i}.maf3.fa.rates >> combined_output.txt
done
shuf -n 28000000 combined_output.txt > combined_output_28M.txt

这两列的具体含义：
第一列：位置的 GERP++ RS（Rejected Substitutions）分数
    这个分数表示该基因组位置的进化约束得分（Constraint Score）。它反映了在比对中此位置受到的 选择压力：
    正分数：表示该位点在进化过程中 非常保守，较少发生取代，表明这个位点可能具有功能重要性。
    负分数：表示该位点相对 不保守，发生了较多的取代，意味着这个位点在进化过程中变化较大，功能性可能不强。
    RS分数的范围：通常是正值（如 0 到 5 或更高），但在某些情况下可以是负值或接近 0。
第二列：物种中实际观察到的取代数量
    这表示在比对序列中实际观察到的 取代数量，通常用来表示某个位置在不同物种中发生了多少次碱基变化。这与 GERP++ 的计算模型中的预期替代次数作对比，帮助估算 RS 分数。
数字越大，表示在进化过程中此位点发生了更多的碱基变化。
数字较小则表示该位点在进化过程中比较保守，变化较少。

library(ggplot2)
data = read.table("combined_output_28M.txt", header=F)
p <- ggplot(data, aes(x = V2)) +
  geom_histogram(binwidth = 0.5, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Distribution of GERP Scores", x = "GERP Score", y = "Frequency") +
  theme_minimal()
pdf("way10_cut_ref_gerp_score_distribution.pdf", height = 6, width = 8)
print(p)
dev.off()

p <- ggplot(data, aes(x = V2)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "Density Plot of GERP Scores", x = "GERP Score", y = "Density") +
  theme_minimal()
pdf("way10_cut_ref_gerp_score_density.pdf", height = 6, width = 8)
print(p)
dev.off()

###计算bed的gerp值
for i in {1..7}
do
    /data/xuebo/software/gerp++KRT/gerpelem -f Chr0${i}.maf3.fa.rates &
done

##统计一下bed文件有多大，覆盖基因组多少区域
for i in {1..7}
do
    awk -v chr="$i" 'BEGIN {OFS="\t"}; {print chr, NR-1, NR, $1, $2}' Chr0${i}.maf3.fa.rates > chr${i}.maf3_clean.fa.bed &
done
## bed file of positive scores
for i in {1..7}
do
    awk '$5 > 0' chr${i}.maf3_clean.fa.bed > chr${i}.maf3_clean.pos.bed
done
##
for i in {1..7}
do
    awk -v chr="Chr$i" 'BEGIN {OFS="\t"} {print chr, $2-1, $3, $4, $5, $6, $7, $8}' Chr0${i}.maf3.fa.rates.elems |sort -k 1,1 -k2,2n   > Chr0${i}.maf3.fa.rates.elems.bed
done
##

```

### Run GERP2
```shell
```bash
/data/xuebo/cucumber/gerp/gerp_cucumber/10_run_gerp
for i in {1..7}
do
    seqtk seq -l0 /data/xuebo/cucumber/gerp/gerp_cucumber/07_combinemaf/msa_fasta/Chr0${i}.maf2.fa | sed 's/> />/g '> Chr0${i}.maf3.fa &
done
for i in {1..7}
do
    /data/xuebo/software/gerp++KRT/gerpcol -f Chr0${i}.maf3.fa -t /data/xuebo/cucumber/gerp/gerp_cucumber/08_neutral_tree/neutral_tree2.txt -v -e ref -j -a &
done
##
for i in {1..7}; do
    awk -v chr="Chr0${i}" '{print chr"\t"$2}' Chr0${i}.maf3.fa.rates >> combined_output.txt
done
##
library(ggplot2)
data = read.table("combined_output.txt", header=F)
p <- ggplot(data, aes(x = V2)) +
  geom_histogram(binwidth = 0.5, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Distribution of GERP Scores", x = "GERP Score", y = "Frequency") +
  theme_minimal()
pdf("gerp_score_distribution.pdf", height = 6, width = 8)
print(p)
dev.off()

p <- ggplot(data, aes(x = V2)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "Density Plot of GERP Scores", x = "GERP Score", y = "Density") +
  theme_minimal()
pdf("gerp_score_density.pdf", height = 6, width = 8)
print(p)
dev.off()

###计算bed的gerp值
for i in {1..7}
do
    /data/xuebo/software/gerp++KRT/gerpelem -f Chr0${i}.maf3.fa.rates &
done

1. region：表示这是一个保守区域的记录，区域的起点和终点由后续列指定。
2. 起始位置（如 1054966）：保守区域在参考基因组中的 起始位置（以基因组坐标为基准，通常从 0 开始）。
3. 终止位置（如 1055487）：保守区域在参考基因组中的 终止位置。
4. 总得分（如 595.633）：保守区域的 RS 分数总和（Rejected Substitutions 总分）。这是该保守区域内所有位点的 RS 分数之和，表示该区域的整体保守性。数值越高，表示区域越保守。
5. P-value（如 2.08002e-07）：表示该保守区域的 显著性 P 值。P 值越低，表明该区域保守性的显著性越高。
6. 元素 ID（如 431054）：这是一个内部 元素 ID，用于唯一标识每个保守区域。这是一个与该区域相关的唯一标识符。
7. 可信度（如 0.5）：保守区域的 可信度（confidence score），通常在 0 到 1 之间。值越接近 1，表示对该区域保守性的置信度越高。
8. 覆盖百分比（如 67.3）：这是保守区域中 保守位点的百分比。表示该区域内多少比例的位点是保守的。百分比越高，说明这个区域内的位点保守性越强。

##统计一下bed文件有多大，覆盖基因组多少区域
for i in {1..7}
do
    awk -v chr="$i" 'BEGIN {OFS="\t"}; {print chr, NR-1, NR, $1, $2}' Chr0${i}.maf3.fa.rates > chr${i}.maf3_clean.fa.bed &
done
## bed file of positive scores
for i in {1..7}
do
    awk '$5 > 0' chr${i}.maf3_clean.fa.bed > chr${i}.maf3_clean.pos.bed
done
##
for i in {1..7}
do
    awk -v chr="Chr$i" 'BEGIN {OFS="\t"} {print chr, $2-1, $3, $4, $5, $6, $7, $8}' Chr0${i}.maf3.fa.rates.elems |sort -k 1,1 -k2,2n   > Chr0${i}.maf3.fa.rates.elems.bed
done
```

### statistical
```shell
```bash
/data/xuebo/cucumber/gerp/gerp_cucumber/11_stat_10way_cut8
cp /data/xuebo/cucumber/annotation/Gene_change_name/AM716.final.gff3 ./
grep -v "chr0" AM716.final.gff3 | grep "CDS" | awk '{sum += ($5 - $4)} END {print sum}'  27155406
grep -v "chr0" AM716.final.gff3 | grep "gene" | awk '{sum += ($5 - $4)} END {print sum}' 77045052
cat ../10_run_gerp/*.maf3.fa.rates.elems.bed > gerp.maf3.fa.rates.elems.bed
awk '{sum += ($3 - $2)} END {print sum}' gerp.maf3.fa.rates.elems.bed 6611055 ##6611055/27155406=0.2434526
sed 's/Chr/chr/g' gerp.maf3.fa.rates.elems.bed | awk '{ if ($2 < 0) $2 = 1; OFS="\t"; print }' > gerp.bed
grep "CDS" AM716.final.gff3 | grep -v "chr0" > ref_CDS.bed
grep "gene" AM716.final.gff3 | grep -v "chr0" > ref_gene.bed
bedtools intersect -a gerp.bed -b ref_CDS.bed > intersect.bed #878671
bedtools intersect -a gerp.bed -b ref_gene.bed > intersect2.bed #2324415

/data/xuebo/cucumber/gerp/gerp_cucumber/11_stat_10way_cutref
cp /data/xuebo/cucumber/annotation/Gene_change_name/AM716.final.gff3 ./
grep -v "chr0" AM716.final.gff3 | grep "CDS" | awk '{sum += ($5 - $4)} END {print sum}'  27155406
grep -v "chr0" AM716.final.gff3 | grep "gene" | awk '{sum += ($5 - $4)} END {print sum}' 77045052
cat ../09_run_gerp/*.maf3.fa.rates.elems.bed > gerp.maf3.fa.rates.elems.bed
awk '{sum += ($3 - $2)} END {print sum}' gerp.maf3.fa.rates.elems.bed 11971616 
sed 's/Chr/chr/g' gerp.maf3.fa.rates.elems.bed > gerp.bed
grep "CDS" AM716.final.gff3 | grep -v "chr0" > ref_CDS.bed
grep "gene" AM716.final.gff3 | grep -v "chr0" > ref_gene.bed
bedtools intersect -a gerp.bed -b ref_CDS.bed > intersect.bed #6532410
bedtools intersect -a gerp.bed -b ref_gene.bed > intersect2.bed #10337740
```

## gerp in watermelon
/data/xuebo/cucumber/gerp/watermelon
### get msa
```shell
```bash
1. zcat ../data/BitterGourd.chr.fa.gz | tail -n +23  | grep -v ">" | sed ':a;N;s/\n/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/g;ta' | sed '1i\>chr0' > chr0_bittergourd.fa
zcat ../data/BitterGourd.chr.fa.gz | head -n 22 | cat - chr0_bittergourd.fa > bittergourd.genome.fa
2. zcat ../data/BottleGourd.chr.fa.gz | tail -n +23  | grep -v ">" | sed ':a;N;s/\n/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/g;ta' | sed '1i\>chr0' > chr0_bottlegourd.fa
zcat ../data/BottleGourd.chr.fa.gz | head -n 22 | sed 's/chr0/chr/g' | cat  - chr0_bottlegourd.fa > bottlegourd.genome.fa
3. zcat ../data/SnakeGourd.chr.fa.gz | tail -n +23  | grep -v ">" | sed ':a;N;s/\n/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/g;ta' | sed '1i\>chr0' > chr0_snakegourd.fa
zcat ../data/SnakeGourd.chr.fa.gz | head -n 22 | sed 's/LG/chr/g' | sed 's/chr0/chr/g' | cat - chr0_snakegourd.fa > snakegourd.genome.fa
4. zcat ../data/Cucumber.chr.fa.gz | tail -n +15  | grep -v ">" | sed ':a;N;s/\n/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/g;ta' | sed '1i\>chr0' > chr0_cucumber.fa
zcat ../data/Cucumber.chr.fa.gz| head -n 14 | cat -  chr0_cucumber.fa > cucumber.genome.fa
5. zcat ../data/Melon.chr.fa.gz | head -n 24 | sed 's/chr0/chr/g' > melon.genome.fa
6. zcat ../data/CA.chr.fa.gz | seqtk seq -l0 | head -n 24 | sed 's/CA02_//g' | sed 's/Chr0/chr/g' | sed 's/Chr/chr/g' > CA.genome.fa
7. zcat ../data/CC.chr.fa.gz | seqtk seq -l0 | head -n 24 | sed 's/CC03_//g' | sed 's/Chr0/chr/g' | sed 's/Chr/chr/g' > CC.genome.fa
8. zcat ../data/CE.chr.fa.gz | seqtk seq -l0 | head -n 24 | sed 's/CE01_//g' | sed 's/Chr0/chr/g' | sed 's/Chr/chr/g' > CE.genome.fa
9. zcat ../data/CLC.chr.fa.gz | seqtk seq -l0 | head -n 24 | sed 's/CLC01_//g' | sed 's/Chr0/chr/g' | sed 's/Chr/chr/g' > CLC.genome.fa
10. zcat ../data/CLV.chr.fa.gz  | seqtk seq -l0 | head -n 24 | sed 's/CLV01_//g' | sed 's/Chr0/chr/g'| sed 's/Chr/chr/g'  > CLV.genome.fa
11. zcat ../data/CM.chr.fa.gz  | seqtk seq -l0 | head -n 24 | sed 's/CM10_//g' | sed 's/Chr0/chr/g' | sed 's/Chr/chr/g' > CM.genome.fa
12. zcat ../data/CN.chr.fa.gz  | seqtk seq -l0 | head -n 24 | sed 's/CN01_//g' | sed 's/Chr0/chr/g' | sed 's/Chr/chr/g' > CN.genome.fa
13. zcat ../data/CR.chr.fa.gz  | seqtk seq -l0 | head -n 24 | sed 's/CR01_//g' | sed 's/Chr0/chr/g' | sed 's/Chr/chr/g' > CR.genome.fa

##Run lastdb.sh on the reference genome 
/data/xuebo/cucumber/gerp/watermelon/ref_CA/00_lastdb
/data/Sunhh/src/align/last/install/v869/bin/lastdb -P 90 -uMAM4 -R01 ./MAM4 /data/xuebo/cucumber/gerp/watermelon/data2/CA.genome.fa &
/data/xuebo/cucumber/gerp/watermelon/ref_CLC/00_lastdb
/data/Sunhh/src/align/last/install/v869/bin/lastdb -P 90 -uMAM4 -R01 ./MAM4 /data/xuebo/cucumber/gerp/watermelon/data2/CLC.genome.fa &
/data/xuebo/cucumber/gerp/watermelon/ref_CLV/00_lastdb
/data/Sunhh/src/align/last/install/v869/bin/lastdb -P 90 -uMAM4 -R01 ./MAM4 /data/xuebo/cucumber/gerp/watermelon/data2/CLV.genome.fa &

##Make alignments
/data/xuebo/cucumber/gerp/watermelon/ref_CA/01_mat
cat ../taxa_12.txt | while read line
do
    /data/Sunhh/src/align/last/install/v869/bin/last-train --revsym --matsym --gapsym -E0.05 -C2 ../00_lastdb/MAM4 /data/xuebo/cucumber/gerp/watermelon/data2/${line}.genome.fa > ref_${line}.mat &
done
/data/xuebo/cucumber/gerp/watermelon/ref_CLC/01_mat
cat ../taxa_12.txt | while read line
do
    /data/Sunhh/src/align/last/install/v869/bin/last-train --revsym --matsym --gapsym -E0.05 -C2 ../00_lastdb/MAM4 /data/xuebo/cucumber/gerp/watermelon/data2/${line}.genome.fa > ref_${line}.mat &
done
/data/xuebo/cucumber/gerp/watermelon/ref_CLV/01_mat
cat ../taxa_12.txt | while read line
do
    /data/Sunhh/src/align/last/install/v869/bin/last-train --revsym --matsym --gapsym -E0.05 -C2 ../00_lastdb/MAM4 /data/xuebo/cucumber/gerp/watermelon/data2/${line}.genome.fa > ref_${line}.mat &
done
###
/data/xuebo/cucumber/gerp/watermelon/ref_CA/02_maf
cat ../taxa_12.txt | while read line
do
    /data/Sunhh/src/align/last/install/v869/bin/lastal -P90 -m50 -E0.05 -C2 -p ../01_mat/ref_${line}.mat ../00_lastdb/MAM4 /data/xuebo/cucumber/gerp/watermelon/data2/${line}.genome.fa > ref_${line}.maf &
done
/data/xuebo/cucumber/gerp/watermelon/ref_CLC/02_maf
cat ../taxa_12.txt | while read line
do
    /data/Sunhh/src/align/last/install/v869/bin/lastal -P90 -m50 -E0.05 -C2 -p ../01_mat/ref_${line}.mat ../00_lastdb/MAM4 /data/xuebo/cucumber/gerp/watermelon/data2/${line}.genome.fa > ref_${line}.maf &
done
/data/xuebo/cucumber/gerp/watermelon/ref_CLC/02_maf
cat ../taxa_12.txt | while read line
do
    /data/Sunhh/src/align/last/install/v869/bin/lastal -P90 -m50 -E0.05 -C2 -p ../01_mat/ref_${line}.mat ../00_lastdb/MAM4 /data/xuebo/cucumber/gerp/watermelon/data2/${line}.genome.fa > ref_${line}.maf &
done
###
/data/xuebo/cucumber/gerp/watermelon/ref_CA/03_axt
cat ../taxa_12.txt | while read line
do
    /data/Sunhh/src/align/last/install/v869/bin/maf-convert axt ../02_maf/ref_${line}.maf > ref_${line}.axt &
done
/data/xuebo/cucumber/gerp/watermelon/ref_CLC/03_axt
cat ../taxa_12.txt | while read line
do
    /data/Sunhh/src/align/last/install/v869/bin/maf-convert axt ../02_maf/ref_${line}.maf > ref_${line}.axt &
done
/data/xuebo/cucumber/gerp/watermelon/ref_CLV/03_axt
cat ../taxa_12.txt | while read line
do
    /data/Sunhh/src/align/last/install/v869/bin/maf-convert axt ../02_maf/ref_${line}.maf > ref_${line}.axt &
done
###
/data/xuebo/cucumber/gerp/watermelon/ref_CA/04_Chain
cat ../taxa_12.txt | while read line
do
    /home/xuebo/bin/x86_64/axtChain ../03_axt/ref_${line}.axt /data/xuebo/cucumber/gerp/watermelon/data2/CA.genome.fa /data/xuebo/cucumber/gerp/watermelon/data2/${line}.genome.fa ref_${line}.chain -linearGap=loose -faQ -faT &
done
nohup /home/xuebo/bin/x86_64/axtChain ../03_axt/ref_CR.axt /data/xuebo/cucumber/gerp/watermelon/data2/CA.genome.fa /data/xuebo/cucumber/gerp/watermelon/data2/CR.genome.fa ref_CR.chain -linearGap=loose -faQ -faT > nohupoutCR 2>& 1 &
nohup /home/xuebo/bin/x86_64/axtChain ../03_axt/ref_snakegourd.axt /data/xuebo/cucumber/gerp/watermelon/data2/CA.genome.fa /data/xuebo/cucumber/gerp/watermelon/data2/snakegourd.genome.fa ref_snakegourd.chain -linearGap=loose -faQ -faT > nohupoutsnakegourd 2>& 1 &
/data/xuebo/cucumber/gerp/watermelon/ref_CLC/04_Chain
cat ../taxa_12.txt | while read line
do
    /home/xuebo/bin/x86_64/axtChain ../03_axt/ref_${line}.axt /data/xuebo/cucumber/gerp/watermelon/data2/CLC.genome.fa /data/xuebo/cucumber/gerp/watermelon/data2/${line}.genome.fa ref_${line}.chain -linearGap=loose -faQ -faT &
done
nohup /home/xuebo/bin/x86_64/axtChain ../03_axt/ref_CR.axt /data/xuebo/cucumber/gerp/watermelon/data2/CLC.genome.fa /data/xuebo/cucumber/gerp/watermelon/data2/CR.genome.fa ref_CR.chain -linearGap=loose -faQ -faT > nohupoutCR 2>& 1 &
nohup /home/xuebo/bin/x86_64/axtChain ../03_axt/ref_snakegourd.axt /data/xuebo/cucumber/gerp/watermelon/data2/CLC.genome.fa /data/xuebo/cucumber/gerp/watermelon/data2/snakegourd.genome.fa ref_snakegourd.chain -linearGap=loose -faQ -faT > nohupoutsnakegourd 2>& 1 &
/data/xuebo/cucumber/gerp/watermelon/ref_CLV/04_Chain
cat ../taxa_12.txt | while read line
do
    /home/xuebo/bin/x86_64/axtChain ../03_axt/ref_${line}.axt /data/xuebo/cucumber/gerp/watermelon/data2/CLV.genome.fa /data/xuebo/cucumber/gerp/watermelon/data2/${line}.genome.fa ref_${line}.chain -linearGap=loose -faQ -faT &
done
nohup /home/xuebo/bin/x86_64/axtChain ../03_axt/ref_snakegourd.axt /data/xuebo/cucumber/gerp/watermelon/data2/CLV.genome.fa /data/xuebo/cucumber/gerp/watermelon/data2/snakegourd.genome.fa ref_snakegourd.chain -linearGap=loose -faQ -faT > nohupoutsnakegourd 2>& 1 &
### 这一步不做了，结果和没有chainMergeSort结果是一模一样的
/data/xuebo/cucumber/gerp/watermelon/ref_CA/04_Chain
cat ../taxa_12.txt | while read line
do
    /home/xuebo/bin/x86_64/chainMergeSort ref_${line}.chain > ref_${line}.merged_chain &
done
```
### CA-getonemaf
```shell
```bash
/data/xuebo/cucumber/gerp/watermelon/ref_CA/05_faSize
/home/xuebo/bin/x86_64/faSize /data/xuebo/cucumber/gerp/watermelon/data2/CA.genome.fa -detailed > ref_size
cat ../taxa_12.txt | while read line
do
    /home/xuebo/bin/x86_64/faSize /data/xuebo/cucumber/gerp/watermelon/data2/${line}.genome.fa -detailed > query_${line}_size &
done
/data/xuebo/cucumber/gerp/watermelon/ref_CA/05_faSize
cat ../taxa_12.txt | while read line
do
    /home/xuebo/bin/x86_64/chainPreNet ../04_Chain/ref_${line}.chain ref_size query_${line}_size ref_${line}.chain_prenet &
done
/data/xuebo/cucumber/gerp/watermelon/ref_CA/05_faSize
cat ../taxa_12.txt | while read line
do
    /home/xuebo/bin/x86_64/chainNet ref_${line}.chain_prenet ref_size query_${line}_size ref_${line}.target_net ref_${line}.query_net &
done
/data/xuebo/cucumber/gerp/watermelon/ref_CA/05_faSize
/home/xuebo/bin/x86_64/faToTwoBit /data/xuebo/cucumber/gerp/watermelon/data2/CA.genome.fa ref_twobit
cat ../taxa_12.txt | while read line
do
    /home/xuebo/bin/x86_64/faToTwoBit /data/xuebo/cucumber/gerp/watermelon/data2/${line}.genome.fa ${line}.query_twobit &
done
/data/xuebo/cucumber/gerp/watermelon/ref_CA/05_faSize
cat ../taxa_12.txt | while read line
do
    /home/xuebo/bin/x86_64/netToAxt ref_${line}.target_net ref_${line}.chain_prenet ref_twobit ${line}.query_twobit ref_${line}.net_axt &
done
/home/xuebo/bin/x86_64/netToAxt ref_snakegourd.target_net ref_snakegourd.chain_prenet ref_twobit snakegourd.query_twobit ref_snakegourd.net_axt & ##这个不能handle，太大了;CR也不能，太大了
/data/xuebo/cucumber/gerp/watermelon/ref_CA/05_faSize
cat ../taxa_12.txt | while read line
do
    /home/xuebo/bin/x86_64/axtToMaf ref_${line}.net_axt ref_size query_${line}_size ref_${line}.net_maf &
done

###
/data/xuebo/cucumber/gerp/watermelon/ref_CA/06_onemaf
cat ../taxa_12.txt | while read line
do
    head -n 29 ../02_maf/ref_${line}.maf > ${line}_net_maf_w_header
    cat ../05_faSize/ref_${line}.net_maf >> ${line}_net_maf_w_header
    /data/Sunhh/src/align/last/install/v869/bin/last-split -m1 ${line}_net_maf_w_header | /data/Sunhh/src/align/last/install/v869/bin/maf-swap | awk -v q="${line}" -v r="ref" '/^s/ {$2 = (++s % 2 ? q "." : r ".") $2} 1' | /data/Sunhh/src/align/last/install/v869/bin/last-split -m1 | /data/Sunhh/src/align/last/install/v869/bin/maf-swap | /data/Sunhh/src/align/last/install/v869/bin/last-postmask > ref_${line}.one_maf &
done

###
rm ref_snakegourd.one_maf
rm ref_CR.one_maf
/data/xuebo/cucumber/gerp/watermelon/ref_CA/07_combinemaf
maf_array=($( ls -d ../06_onemaf/*.one_maf ))
combined_maf=./combined.maf
cat ${maf_array[@]:0:1} | sed -n '/##maf version=1 scoring=blastz/,$p' > \
${maf_array[@]:0:1}_tmp
cat ${maf_array[@]:1:1} | sed -n '/##maf version=1 scoring=blastz/,$p' > \
${maf_array[@]:1:1}_tmp
/data/xuebo/software/multiz-20190527/multiz ${maf_array[@]:0:1}_tmp ${maf_array[@]:1:1}_tmp 1 > $combined_maf
for maf in ${maf_array[@]:2};
do
  echo "processing " $maf
 cat $maf | sed -n '/##maf version=1 scoring=blastz/,$p' > \
 "$maf"_tmp
 /data/xuebo/software/multiz-20190527/multiz $combined_maf "$maf"_tmp 1 > "$combined_maf"_tmp
 mv "$combined_maf"_tmp $combined_maf
done
# and filter mafs so all blocks have Zea mays and are at least 20 bp long
/home/xuebo/bin/x86_64/mafFilter -minCol=20 -needComp="ref" combined.maf > combined.maf.filtered ##这个是用来计算gerp值的
/home/xuebo/bin/x86_64/mafFilter -minCol=20 -minRow=11 -needComp="ref" combined.maf > combined.maf2.filtered ##这个事用来计算tree的

###
/data/xuebo/cucumber/gerp/watermelon/ref_CA/07_combinemaf
##Convert maf files to fasta files
mkdir split_maf
outroot=./split_maf/
/home/xuebo/bin/x86_64/mafSplit -byTarget dummy.bed $outroot combined.maf.filtered -useFullSequenceName
mkdir ref_split
/home/xuebo/bin/x86_64/faSplit byname /data/xuebo/cucumber/gerp/watermelon/data2/CA.genome.fa ./ref_split/
for i in {0..11}
do
    mv chr${i}.fa chr${i}.maf.fa
done
##
mkdir msa_fasta
maf_dir=($( ls -d /data/xuebo/cucumber/gerp/watermelon/ref_CA/07_combinemaf/split_maf/*chr*.maf ))
for maf_file in "${maf_dir[@]}"
do
  chr=$(basename $maf_file | sed -e 's/\.fa$//')

  fasta=./msa_fasta/"$chr".fa
  ref_rm_chr=./ref_split/"$chr".fa
  rm_fasta=./msa_fasta/"$chr"_rm.fa
  /data/zhaojiantao/tools/phast/bin/msa_view $maf_file -f -G 1 --refseq $ref_rm_chr --missing-as-indels > $fasta
  sed -E 's/> />/g' $fasta > "$fasta"_tmp && mv "$fasta"_tmp $fasta
  perl /data/xuebo/cucumber/gerp/From_HuffordLab/matchMasking.pl \
  --ref $ref_rm_chr \
  --fasta $fasta \
  --out $rm_fasta
done
```

### CA-tree
```shell
```bash
cp ../07_combinemaf/combined.maf2.filtered ./
mkdir split_maf
outroot=./split_maf/
/home/xuebo/bin/x86_64/mafSplit -byTarget dummy.bed $outroot combined.maf2.filtered -useFullSequenceName
mkdir ref_split
/home/xuebo/bin/x86_64/faSplit byname /data/xuebo/cucumber/gerp/watermelon/data2/CA.genome.fa ./ref_split/
for i in {0..11}
do
    mv chr${i}.fa chr${i}.maf.fa
done
##
mkdir msa_fasta
maf_dir=($( ls -d /data/xuebo/cucumber/gerp/watermelon/ref_CA/08_tree/split_maf/*chr*.maf ))
for maf_file in "${maf_dir[@]}"
do
  chr=$(basename $maf_file | sed -e 's/\.fa$//')

  fasta=./msa_fasta/"$chr".fa
  ref_rm_chr=./ref_split/"$chr".fa
  rm_fasta=./msa_fasta/"$chr"_rm.fa
  /data/zhaojiantao/tools/phast/bin/msa_view $maf_file -f -G 1 --refseq $ref_rm_chr --missing-as-indels > $fasta
  sed -E 's/> />/g' $fasta > "$fasta"_tmp && mv "$fasta"_tmp $fasta
  perl /data/xuebo/cucumber/gerp/From_HuffordLab/matchMasking.pl \
  --ref $ref_rm_chr \
  --fasta $fasta \
  --out $rm_fasta
done
for i in {1..11}
do
    /data/zhaojiantao/tools/phast/bin/msa_view chr${i}.maf.fa -f --gap-strip 11  > chr${i}.maf2.fa &
done
for i in {1..11}
do
    seqtk seq -l0 chr${i}.maf2.fa > chr${i}.maf3.fa &
done
###合并染色体
#!/bin/bash
# 定义包含所有染色体fasta文件的列表
fasta_files=("chr1.maf3.fa" "chr2.maf3.fa" "chr3.maf3.fa" "chr4.maf3.fa" "chr5.maf3.fa" "chr6.maf3.fa"  "chr7.maf3.fa" "chr8.maf3.fa" "chr9.maf3.fa" "chr10.maf3.fa" "chr11.maf3.fa")  # 修改为您的fasta文件名
# 定义输出文件
output_file="merged_species.fasta"
# 清空或创建输出文件
> $output_file
# 使用一个关联数组存储每个物种的合并序列
declare -A species_sequences
# 循环处理每个fasta文件
for fasta_file in "${fasta_files[@]}"; do
    # 读取每个fasta文件中的内容
    while read line; do
        # 如果是物种行（>开头的行）
        if [[ $line == ">"* ]]; then
            # 提取物种名称（去掉>）
            species=$(echo "$line" | sed 's/>//')
        else
            # 将序列添加到对应的物种
            species_sequences[$species]+="$line"
        fi
    done < "$fasta_file"
done
# 将合并后的结果写入到输出文件
for species in "${!species_sequences[@]}"; do
    echo ">$species" >> $output_file
    echo "${species_sequences[$species]}" >> $output_file
done
echo "合并后的序列已写入 $output_file"

###随机挑选10k
# vim select_random_msa.py
import random
from Bio import AlignIO
# 输入的 MSA 文件（fasta 格式）
msa_file = "merged_species.fasta"  # 替换为您的MSA文件名
output_file = "msa_random_10k.fasta"
num_bases = 10000  # 需要随机选择的碱基数量
# 读取MSA文件
alignment = AlignIO.read(msa_file, "fasta")
# 获取对齐的总长度
alignment_length = alignment.get_alignment_length()
# 确保 alignment_length 足够长，可以抽取10,000个碱基
if alignment_length < num_bases:
    raise ValueError(f"对齐长度为 {alignment_length} bp，无法挑选 {num_bases} 个碱基")
# 随机选择 num_bases 个位置
random_positions = sorted(random.sample(range(alignment_length), num_bases))
# 创建新的序列容器
selected_sequences = {}
# 遍历每个物种，并从随机选择的位置提取碱基
for record in alignment:
    sequence = ""
    for pos in random_positions:
        sequence += record.seq[pos]
    selected_sequences[record.id] = sequence
# 将随机选择的子序列写入输出fasta文件
with open(output_file, "w") as out_fasta:
    for species, sequence in selected_sequences.items():
        out_fasta.write(f">{species}\n")
        out_fasta.write(f"{sequence}\n")
print(f"随机选择的 {num_bases} 个碱基已写入 {output_file}")
```
### CA-gerp
```shell
```bash
for i in {0..11}
do
    seqtk seq -l0 /data/xuebo/cucumber/gerp/watermelon/ref_CA/07_combinemaf/msa_fasta/chr${i}.maf.fa > chr${i}.maf3.fa &
done
for i in {0..11}
do
    /data/xuebo/software/gerp++KRT/gerpcol -f chr${i}.maf3.fa -t /data/xuebo/cucumber/gerp/watermelon/ref_CA/08_tree/neutral_tree.txt -v -e ref -j -a &
done
## /data/zhaojiantao/tools/phast/bin/msa_view chr0.maf.fa -f --gap-strip 1 > xx 
##
for i in {0..11}; do
    awk -v chr="chr${i}" '{print chr"\t"$2}' chr${i}.maf3.fa.rates >> combined_output.txt
done
shuf -n 3000000 combined_output.txt > combined_output_shuf.txt
##
library(ggplot2)
data = read.table("combined_output_shuf.txt", header=F)
p <- ggplot(data, aes(x = V2)) +
  geom_histogram(binwidth = 0.05, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "refCA: Distribution of GERP Scores", x = "GERP Score", y = "Frequency") +
  theme_minimal()
pdf("refCA_gerp_score_distribution.pdf", height = 6, width = 8)
print(p)
dev.off()

p <- ggplot(data, aes(x = V2)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "refCA: Density Plot of GERP Scores", x = "GERP Score", y = "Density") +
  theme_minimal()
pdf("refCA_gerp_score_density.pdf", height = 6, width = 8)
print(p)
dev.off()
```

### CLC-getonemaf
```shell
```bash
/data/xuebo/cucumber/gerp/watermelon/ref_CLC/05_faSize
/home/xuebo/bin/x86_64/faSize /data/xuebo/cucumber/gerp/watermelon/data2/CLC.genome.fa -detailed > ref_size
cat ../taxa_12.txt | while read line
do
    /home/xuebo/bin/x86_64/faSize /data/xuebo/cucumber/gerp/watermelon/data2/${line}.genome.fa -detailed > query_${line}_size &
done
/data/xuebo/cucumber/gerp/watermelon/ref_CLC/05_faSize
cat ../taxa_12.txt | while read line
do
    /home/xuebo/bin/x86_64/chainPreNet ../04_Chain/ref_${line}.chain ref_size query_${line}_size ref_${line}.chain_prenet &
done
/data/xuebo/cucumber/gerp/watermelon/ref_CLC/05_faSize
cat ../taxa_12.txt | while read line
do
    /home/xuebo/bin/x86_64/chainNet ref_${line}.chain_prenet ref_size query_${line}_size ref_${line}.target_net ref_${line}.query_net &
done
/data/xuebo/cucumber/gerp/watermelon/ref_CLC/05_faSize
/home/xuebo/bin/x86_64/faToTwoBit /data/xuebo/cucumber/gerp/watermelon/data2/CLC.genome.fa ref_twobit
cat ../taxa_12.txt | while read line
do
    /home/xuebo/bin/x86_64/faToTwoBit /data/xuebo/cucumber/gerp/watermelon/data2/${line}.genome.fa ${line}.query_twobit &
done
/data/xuebo/cucumber/gerp/watermelon/ref_CLC/05_faSize
cat ../taxa_12.txt | while read line
do
    /home/xuebo/bin/x86_64/netToAxt ref_${line}.target_net ref_${line}.chain_prenet ref_twobit ${line}.query_twobit ref_${line}.net_axt &
done
/home/xuebo/bin/x86_64/netToAxt ref_snakegourd.target_net ref_snakegourd.chain_prenet ref_twobit snakegourd.query_twobit ref_snakegourd.net_axt &
/data/xuebo/cucumber/gerp/watermelon/ref_CLC/05_faSize
cat ../taxa_12.txt | while read line
do
    /home/xuebo/bin/x86_64/axtToMaf ref_${line}.net_axt ref_size query_${line}_size ref_${line}.net_maf &
done

###
/data/xuebo/cucumber/gerp/watermelon/ref_CLC/06_onemaf
cat ../taxa_12.txt | while read line
do
    head -n 29 ../02_maf/ref_${line}.maf > ${line}_net_maf_w_header
    cat ../05_faSize/ref_${line}.net_maf >> ${line}_net_maf_w_header
    /data/Sunhh/src/align/last/install/v869/bin/last-split -m1 ${line}_net_maf_w_header | /data/Sunhh/src/align/last/install/v869/bin/maf-swap | awk -v q="${line}" -v r="ref" '/^s/ {$2 = (++s % 2 ? q "." : r ".") $2} 1' | /data/Sunhh/src/align/last/install/v869/bin/last-split -m1 | /data/Sunhh/src/align/last/install/v869/bin/maf-swap | /data/Sunhh/src/align/last/install/v869/bin/last-postmask > ref_${line}.one_maf &
done

###
rm ref_snakegourd.one_maf
/data/xuebo/cucumber/gerp/watermelon/ref_CLC/07_combinemaf
maf_array=($( ls -d ../06_onemaf/*.one_maf ))
combined_maf=./combined.maf
cat ${maf_array[@]:0:1} | sed -n '/##maf version=1 scoring=blastz/,$p' > \
${maf_array[@]:0:1}_tmp
cat ${maf_array[@]:1:1} | sed -n '/##maf version=1 scoring=blastz/,$p' > \
${maf_array[@]:1:1}_tmp
/data/xuebo/software/multiz-20190527/multiz ${maf_array[@]:0:1}_tmp ${maf_array[@]:1:1}_tmp 1 > $combined_maf
for maf in ${maf_array[@]:2};
do
  echo "processing " $maf
 cat $maf | sed -n '/##maf version=1 scoring=blastz/,$p' > \
 "$maf"_tmp
 /data/xuebo/software/multiz-20190527/multiz $combined_maf "$maf"_tmp 1 > "$combined_maf"_tmp
 mv "$combined_maf"_tmp $combined_maf
done
# and filter mafs so all blocks have Zea mays and are at least 20 bp long
/home/xuebo/bin/x86_64/mafFilter -minCol=20 -needComp="ref" combined.maf > combined.maf.filtered ##这个是用来计算gerp值的
/home/xuebo/bin/x86_64/mafFilter -minCol=20 -minRow=11 -needComp="ref" combined.maf > combined.maf2.filtered ##这个事用来计算tree的

###
/data/xuebo/cucumber/gerp/watermelon/ref_CLC/07_combinemaf
##Convert maf files to fasta files
mkdir split_maf
outroot=./split_maf/
/home/xuebo/bin/x86_64/mafSplit -byTarget dummy.bed $outroot combined.maf.filtered -useFullSequenceName
mkdir ref_split
/home/xuebo/bin/x86_64/faSplit byname /data/xuebo/cucumber/gerp/watermelon/data2/CLC.genome.fa ./ref_split/
for i in {0..11}
do
    mv chr${i}.fa chr${i}.maf.fa
done
##
mkdir msa_fasta
maf_dir=($( ls -d /data/xuebo/cucumber/gerp/watermelon/ref_CLC/07_combinemaf/split_maf/*chr*.maf ))
for maf_file in "${maf_dir[@]}"
do
  chr=$(basename $maf_file | sed -e 's/\.fa$//')

  fasta=./msa_fasta/"$chr".fa
  ref_rm_chr=./ref_split/"$chr".fa
  rm_fasta=./msa_fasta/"$chr"_rm.fa
  /data/zhaojiantao/tools/phast/bin/msa_view $maf_file -f -G 1 --refseq $ref_rm_chr  > $fasta
  sed -E 's/> />/g' $fasta > "$fasta"_tmp && mv "$fasta"_tmp $fasta
  perl /data/xuebo/cucumber/gerp/From_HuffordLab/matchMasking.pl \
  --ref $ref_rm_chr \
  --fasta $fasta \
  --out $rm_fasta
done
```

### CLC-tree
```shell
```bash
cp ../07_combinemaf/combined.maf2.filtered ./
mkdir split_maf
outroot=./split_maf/
/home/xuebo/bin/x86_64/mafSplit -byTarget dummy.bed $outroot combined.maf2.filtered -useFullSequenceName
mkdir ref_split
/home/xuebo/bin/x86_64/faSplit byname /data/xuebo/cucumber/gerp/watermelon/data2/CLC.genome.fa ./ref_split/
for i in {0..11}
do
    mv chr${i}.fa chr${i}.maf.fa
done
##
mkdir msa_fasta
maf_dir=($( ls -d /data/xuebo/cucumber/gerp/watermelon/ref_CLC/08_tree/split_maf/*chr*.maf ))
for maf_file in "${maf_dir[@]}"
do
  chr=$(basename $maf_file | sed -e 's/\.fa$//')

  fasta=./msa_fasta/"$chr".fa
  ref_rm_chr=./ref_split/"$chr".fa
  rm_fasta=./msa_fasta/"$chr"_rm.fa
  /data/zhaojiantao/tools/phast/bin/msa_view $maf_file -f -G 1 --refseq $ref_rm_chr --missing-as-indels > $fasta
  sed -E 's/> />/g' $fasta > "$fasta"_tmp && mv "$fasta"_tmp $fasta
  perl /data/xuebo/cucumber/gerp/From_HuffordLab/matchMasking.pl \
  --ref $ref_rm_chr \
  --fasta $fasta \
  --out $rm_fasta
done
for i in {1..11}
do
    /data/zhaojiantao/tools/phast/bin/msa_view chr${i}.maf.fa -f --gap-strip 11  > chr${i}.maf2.fa &
done
for i in {1..11}
do
    seqtk seq -l0 chr${i}.maf2.fa > chr${i}.maf3.fa &
done
###合并染色体
#!/bin/bash
# 定义包含所有染色体fasta文件的列表
fasta_files=("chr1.maf3.fa" "chr2.maf3.fa" "chr3.maf3.fa" "chr4.maf3.fa" "chr5.maf3.fa" "chr6.maf3.fa"  "chr7.maf3.fa" "chr8.maf3.fa" "chr9.maf3.fa" "chr10.maf3.fa" "chr11.maf3.fa")  # 修改为您的fasta文件名
# 定义输出文件
output_file="merged_species.fasta"
# 清空或创建输出文件
> $output_file
# 使用一个关联数组存储每个物种的合并序列
declare -A species_sequences
# 循环处理每个fasta文件
for fasta_file in "${fasta_files[@]}"; do
    # 读取每个fasta文件中的内容
    species=""
    while read line; do
        # 如果是物种行（>开头的行）
        if [[ $line == ">"* ]]; then
            # 提取物种名称（去掉>），并重置species
            species=$(echo "$line" | sed 's/>//')
        elif [[ -n "$species" ]]; then  # 只有species已经被设置时才处理序列行
            # 将序列添加到对应的物种
            species_sequences[$species]+="$line"
        fi
    done < "$fasta_file"
done
# 将合并后的结果写入到输出文件
for species in "${!species_sequences[@]}"; do
    echo ">$species" >> $output_file
    echo "${species_sequences[$species]}" >> $output_file
done
echo "合并后的序列已写入 $output_file"


###随机挑选10k
# vim select_random_msa.py
import random
from Bio import AlignIO
# 输入的 MSA 文件（fasta 格式）
msa_file = "merged_species.fasta"  # 替换为您的MSA文件名
output_file = "msa_random_10k.fasta"
num_bases = 10000  # 需要随机选择的碱基数量
# 读取MSA文件
alignment = AlignIO.read(msa_file, "fasta")
# 获取对齐的总长度
alignment_length = alignment.get_alignment_length()
# 确保 alignment_length 足够长，可以抽取10,000个碱基
if alignment_length < num_bases:
    raise ValueError(f"对齐长度为 {alignment_length} bp，无法挑选 {num_bases} 个碱基")
# 随机选择 num_bases 个位置
random_positions = sorted(random.sample(range(alignment_length), num_bases))
# 创建新的序列容器
selected_sequences = {}
# 遍历每个物种，并从随机选择的位置提取碱基
for record in alignment:
    sequence = ""
    for pos in random_positions:
        sequence += record.seq[pos]
    selected_sequences[record.id] = sequence
# 将随机选择的子序列写入输出fasta文件
with open(output_file, "w") as out_fasta:
    for species, sequence in selected_sequences.items():
        out_fasta.write(f">{species}\n")
        out_fasta.write(f"{sequence}\n")
print(f"随机选择的 {num_bases} 个碱基已写入 {output_file}")
```

### CLC-gerp
```shell
```bash
for i in {0..11}
do
    seqtk seq -l0 /data/xuebo/cucumber/gerp/watermelon/ref_CLC/07_combinemaf/msa_fasta/chr${i}.maf.fa |  sed 's/*/-/g' > chr${i}.maf3.fa &
done
head -n 2 chr0.maf3.fa | seqtk comp 
head -n 2 chr0.maf3.fa | tail +2 | fold -w1 | sort | uniq -c
seqtk seq -l0 /data/xuebo/cucumber/gerp/watermelon/ref_CLC/07_combinemaf/msa_fasta/chr0.maf.fa | head -n 2 | seqtk comp 
seqtk seq -l0 /data/xuebo/cucumber/gerp/watermelon/ref_CLC/07_combinemaf/msa_fasta/chr0.maf.fa  | head -n 2 | tail +2 | fold -w1 | sort | uniq -c
sed '/^>/! s/n/A/g' chr0.maf3.fa > xx
/data/xuebo/software/gerp++KRT/gerpcol -f chr0.maf3.fa -t /data/xuebo/cucumber/gerp/watermelon/ref_CLC/08_tree/neutral_tree.txt -v -e ref -j -a
/data/xuebo/software/gerp++KRT/gerpcol -f xx -t /data/xuebo/cucumber/gerp/watermelon/ref_CLC/08_tree/neutral_tree.txt -v -e ref -j -a
######由此可以看出，ref是N的，结果是算不出来的
for i in {0..11}
do
    /data/xuebo/software/gerp++KRT/gerpcol -f chr${i}.maf3.fa -t /data/xuebo/cucumber/gerp/watermelon/ref_CLC/08_tree/neutral_tree.txt -v -e ref -j -a &
done
## /data/zhaojiantao/tools/phast/bin/msa_view chr0.maf.fa -f --gap-strip 1 > xx 
##
for i in {0..11}; do
    awk -v chr="chr${i}" '{print chr"\t"$2}' chr${i}.maf3.fa.rates >> combined_output.txt
done
shuf -n 3000000 combined_output.txt > combined_output_shuf.txt
##
library(ggplot2)
data = read.table("combined_output_shuf.txt", header=F)
p <- ggplot(data, aes(x = V2)) +
  geom_histogram(binwidth = 0.05, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "refCLC: Distribution of GERP Scores", x = "GERP Score", y = "Frequency") +
  theme_minimal()
pdf("refCLC_gerp_score_distribution.pdf", height = 6, width = 8)
print(p)
dev.off()

p <- ggplot(data, aes(x = V2)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "refCLC: Density Plot of GERP Scores", x = "GERP Score", y = "Density") +
  theme_minimal()
pdf("refCLC_gerp_score_density.pdf", height = 6, width = 8)
print(p)
dev.off()
```

### CLV-getonemaf
```shell
```bash
/data/xuebo/cucumber/gerp/watermelon/ref_CLV/05_faSize
/home/xuebo/bin/x86_64/faSize /data/xuebo/cucumber/gerp/watermelon/data2/CLV.genome.fa -detailed > ref_size
cat ../taxa_12.txt | while read line
do
    /home/xuebo/bin/x86_64/faSize /data/xuebo/cucumber/gerp/watermelon/data2/${line}.genome.fa -detailed > query_${line}_size &
done
/data/xuebo/cucumber/gerp/watermelon/ref_CLV/05_faSize
cat ../taxa_12.txt | while read line
do
    /home/xuebo/bin/x86_64/chainPreNet ../04_Chain/ref_${line}.chain ref_size query_${line}_size ref_${line}.chain_prenet &
done
/data/xuebo/cucumber/gerp/watermelon/ref_CLV/05_faSize
cat ../taxa_12.txt | while read line
do
    /home/xuebo/bin/x86_64/chainNet ref_${line}.chain_prenet ref_size query_${line}_size ref_${line}.target_net ref_${line}.query_net &
done                 
/data/xuebo/cucumber/gerp/watermelon/ref_CLV/05_faSize
/home/xuebo/bin/x86_64/faToTwoBit /data/xuebo/cucumber/gerp/watermelon/data2/CLV.genome.fa ref_twobit
cat ../taxa_12.txt | while read line
do
    /home/xuebo/bin/x86_64/faToTwoBit /data/xuebo/cucumber/gerp/watermelon/data2/${line}.genome.fa ${line}.query_twobit &
done
/data/xuebo/cucumber/gerp/watermelon/ref_CLV/05_faSize
cat ../taxa_12.txt | while read line
do
    /home/xuebo/bin/x86_64/netToAxt ref_${line}.target_net ref_${line}.chain_prenet ref_twobit ${line}.query_twobit ref_${line}.net_axt &
done
/home/xuebo/bin/x86_64/netToAxt ref_snakegourd.target_net ref_snakegourd.chain_prenet ref_twobit snakegourd.query_twobit ref_snakegourd.net_axt & ##这个不能handle，太大了
/data/xuebo/cucumber/gerp/watermelon/ref_CLV/05_faSize
cat ../taxa_12.txt | while read line
do
    /home/xuebo/bin/x86_64/axtToMaf ref_${line}.net_axt ref_size query_${line}_size ref_${line}.net_maf &
done

###
/data/xuebo/cucumber/gerp/watermelon/ref_CLV/06_onemaf
cat ../taxa_12.txt | while read line
do
    head -n 29 ../02_maf/ref_${line}.maf > ${line}_net_maf_w_header
    cat ../05_faSize/ref_${line}.net_maf >> ${line}_net_maf_w_header
    /data/Sunhh/src/align/last/install/v869/bin/last-split -m1 ${line}_net_maf_w_header | /data/Sunhh/src/align/last/install/v869/bin/maf-swap | awk -v q="${line}" -v r="ref" '/^s/ {$2 = (++s % 2 ? q "." : r ".") $2} 1' | /data/Sunhh/src/align/last/install/v869/bin/last-split -m1 | /data/Sunhh/src/align/last/install/v869/bin/maf-swap | /data/Sunhh/src/align/last/install/v869/bin/last-postmask > ref_${line}.one_maf &
done

###
rm ref_snakegourd.one_maf
/data/xuebo/cucumber/gerp/watermelon/ref_CLV/07_combinemaf
maf_array=($( ls -d ../06_onemaf/*.one_maf ))
combined_maf=./combined.maf
cat ${maf_array[@]:0:1} | sed -n '/##maf version=1 scoring=blastz/,$p' > \
${maf_array[@]:0:1}_tmp
cat ${maf_array[@]:1:1} | sed -n '/##maf version=1 scoring=blastz/,$p' > \
${maf_array[@]:1:1}_tmp
/data/xuebo/software/multiz-20190527/multiz ${maf_array[@]:0:1}_tmp ${maf_array[@]:1:1}_tmp 1 > $combined_maf
for maf in ${maf_array[@]:2};
do
  echo "processing " $maf
 cat $maf | sed -n '/##maf version=1 scoring=blastz/,$p' > \
 "$maf"_tmp
 /data/xuebo/software/multiz-20190527/multiz $combined_maf "$maf"_tmp 1 > "$combined_maf"_tmp
 mv "$combined_maf"_tmp $combined_maf
done
# and filter mafs so all blocks have Zea mays and are at least 20 bp long
/home/xuebo/bin/x86_64/mafFilter -minCol=20 -needComp="ref" combined.maf > combined.maf.filtered ##这个是用来计算gerp值的
/home/xuebo/bin/x86_64/mafFilter -minCol=20 -minRow=11 -needComp="ref" combined.maf > combined.maf2.filtered ##这个事用来计算tree的

###
/data/xuebo/cucumber/gerp/watermelon/ref_CA/07_combinemaf
##Convert maf files to fasta files
mkdir split_maf
outroot=./split_maf/
/home/xuebo/bin/x86_64/mafSplit -byTarget dummy.bed $outroot combined.maf.filtered -useFullSequenceName
mkdir ref_split
/home/xuebo/bin/x86_64/faSplit byname /data/xuebo/cucumber/gerp/watermelon/data2/CLV.genome.fa ./ref_split/
for i in {0..11}
do
    mv chr${i}.fa chr${i}.maf.fa
done
##
mkdir msa_fasta
maf_dir=($( ls -d /data/xuebo/cucumber/gerp/watermelon/ref_CLV/07_combinemaf/split_maf/*chr*.maf ))
for maf_file in "${maf_dir[@]}"
do
  chr=$(basename $maf_file | sed -e 's/\.fa$//')

  fasta=./msa_fasta/"$chr".fa
  ref_rm_chr=./ref_split/"$chr".fa
  rm_fasta=./msa_fasta/"$chr"_rm.fa
  /data/zhaojiantao/tools/phast/bin/msa_view $maf_file -f -G 1 --refseq $ref_rm_chr --missing-as-indels > $fasta
  sed -E 's/> />/g' $fasta > "$fasta"_tmp && mv "$fasta"_tmp $fasta
  perl /data/xuebo/cucumber/gerp/From_HuffordLab/matchMasking.pl \
  --ref $ref_rm_chr \
  --fasta $fasta \
  --out $rm_fasta
done
```

### CLV-tree
```shell
```bash
cp ../07_combinemaf/combined.maf2.filtered ./
mkdir split_maf
outroot=./split_maf/
/home/xuebo/bin/x86_64/mafSplit -byTarget dummy.bed $outroot combined.maf2.filtered -useFullSequenceName
mkdir ref_split
/home/xuebo/bin/x86_64/faSplit byname /data/xuebo/cucumber/gerp/watermelon/data2/CLV.genome.fa ./ref_split/
for i in {0..11}
do
    mv chr${i}.fa chr${i}.maf.fa
done
##
mkdir msa_fasta
maf_dir=($( ls -d /data/xuebo/cucumber/gerp/watermelon/ref_CLV/08_tree/split_maf/*chr*.maf ))
for maf_file in "${maf_dir[@]}"
do
  chr=$(basename $maf_file | sed -e 's/\.fa$//')

  fasta=./msa_fasta/"$chr".fa
  ref_rm_chr=./ref_split/"$chr".fa
  rm_fasta=./msa_fasta/"$chr"_rm.fa
  /data/zhaojiantao/tools/phast/bin/msa_view $maf_file -f -G 1 --refseq $ref_rm_chr --missing-as-indels > $fasta
  sed -E 's/> />/g' $fasta > "$fasta"_tmp && mv "$fasta"_tmp $fasta
  perl /data/xuebo/cucumber/gerp/From_HuffordLab/matchMasking.pl \
  --ref $ref_rm_chr \
  --fasta $fasta \
  --out $rm_fasta
done
for i in {1..11}
do
    /data/zhaojiantao/tools/phast/bin/msa_view chr${i}.maf.fa -f --gap-strip 11  > chr${i}.maf2.fa &
done
for i in {1..11}
do
    seqtk seq -l0 chr${i}.maf2.fa > chr${i}.maf3.fa &
done
###合并染色体
#!/bin/bash
# 定义包含所有染色体fasta文件的列表
fasta_files=("chr1.maf3.fa" "chr2.maf3.fa" "chr3.maf3.fa" "chr4.maf3.fa" "chr5.maf3.fa" "chr6.maf3.fa"  "chr7.maf3.fa" "chr8.maf3.fa" "chr9.maf3.fa" "chr10.maf3.fa" "chr11.maf3.fa")  # 修改为您的fasta文件名
# 定义输出文件
output_file="merged_species.fasta"
# 清空或创建输出文件
> $output_file
# 使用一个关联数组存储每个物种的合并序列
declare -A species_sequences
# 循环处理每个fasta文件
for fasta_file in "${fasta_files[@]}"; do
    # 读取每个fasta文件中的内容
    while read line; do
        # 如果是物种行（>开头的行）
        if [[ $line == ">"* ]]; then
            # 提取物种名称（去掉>）
            species=$(echo "$line" | sed 's/>//')
        else
            # 将序列添加到对应的物种
            species_sequences[$species]+="$line"
        fi
    done < "$fasta_file"
done
# 将合并后的结果写入到输出文件
for species in "${!species_sequences[@]}"; do
    echo ">$species" >> $output_file
    echo "${species_sequences[$species]}" >> $output_file
done
echo "合并后的序列已写入 $output_file"

###随机挑选10k
# vim select_random_msa.py
import random
from Bio import AlignIO
# 输入的 MSA 文件（fasta 格式）
msa_file = "merged_species.fasta"  # 替换为您的MSA文件名
output_file = "msa_random_10k.fasta"
num_bases = 10000  # 需要随机选择的碱基数量
# 读取MSA文件
alignment = AlignIO.read(msa_file, "fasta")
# 获取对齐的总长度
alignment_length = alignment.get_alignment_length()
# 确保 alignment_length 足够长，可以抽取10,000个碱基
if alignment_length < num_bases:
    raise ValueError(f"对齐长度为 {alignment_length} bp，无法挑选 {num_bases} 个碱基")
# 随机选择 num_bases 个位置
random_positions = sorted(random.sample(range(alignment_length), num_bases))
# 创建新的序列容器
selected_sequences = {}
# 遍历每个物种，并从随机选择的位置提取碱基
for record in alignment:
    sequence = ""
    for pos in random_positions:
        sequence += record.seq[pos]
    selected_sequences[record.id] = sequence
# 将随机选择的子序列写入输出fasta文件
with open(output_file, "w") as out_fasta:
    for species, sequence in selected_sequences.items():
        out_fasta.write(f">{species}\n")
        out_fasta.write(f"{sequence}\n")
print(f"随机选择的 {num_bases} 个碱基已写入 {output_file}")
```
### CLV-gerp
```shell
```bash
for i in {0..11}
do
    seqtk seq -l0 /data/xuebo/cucumber/gerp/watermelon/ref_CLV/07_combinemaf/msa_fasta/chr${i}.maf.fa > chr${i}.maf3.fa &
done
for i in {0..11}
do
    /data/xuebo/software/gerp++KRT/gerpcol -f chr${i}.maf3.fa -t /data/xuebo/cucumber/gerp/watermelon/ref_CLV/08_tree/neutral_tree.txt -v -e ref -j -a &
done
## /data/zhaojiantao/tools/phast/bin/msa_view chr0.maf.fa -f --gap-strip 1 > xx 
for i in {0..11}; do
    awk -v chr="chr${i}" '{print chr"\t"$2}' chr${i}.maf3.fa.rates >> combined_output.txt
done
shuf -n 3000000 combined_output.txt > combined_output_shuf.txt
##
library(ggplot2)
data = read.table("combined_output_shuf.txt", header=F)
p <- ggplot(data, aes(x = V2)) +
  geom_histogram(binwidth = 0.05, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "refCLV: Distribution of GERP Scores", x = "GERP Score", y = "Frequency") +
  theme_minimal()
pdf("refCLV_gerp_score_distribution.pdf", height = 6, width = 8)
print(p)
dev.off()

p <- ggplot(data, aes(x = V2)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "refCLV: Density Plot of GERP Scores", x = "GERP Score", y = "Density") +
  theme_minimal()
pdf("refCLV_gerp_score_density.pdf", height = 6, width = 8)
print(p)
dev.off()
```

### 统计比对质量
```shell
```bash
bedtools makewindows -g genomesize.txt -w 10000 > intervals_10k.bed
awk '{print $1"\t"$2}' SNP_final.geno  | sort -k1,1V -k2,2n | awk '{print $1"\t"$2"\t"$2+1}' > SNP_final.bed
bedtools map -a intervals_10k.bed -b SNP_final.bed -c 2 -o count > intervals_with_count.bed
bedtools makewindows -g genomesize.txt -w 100000 > intervals_100k.bed
bedtools map -a intervals_100k.bed -b SNP_final.bed -c 2 -o count > intervals_with_count_100k.bed
#######算一下比对的质量
grep "^s" chr3.maf | awk '{print $2, $3, $4, $5}' | grep "ref" | less -S
grep "score=" chr3.maf | cut -d"=" -f2 |  cut -d"." -f1 > xx
#!/bin/bash
# 输出表头
echo -e "Chromosome\tStart\tScore\tLength\tSpecies_Count"
# 使用awk解析MAF文件
awk '
BEGIN {
    FS = "[ ]+"
    OFS = "\t"
}
# 捕获a score行，并存储比对块的评分
/^a score/ {
    # 如果已经有了一个比对块，则输出它
    if (species_count > 0) {
        print chr, start, score, len, species_count
    }
    # 处理新的比对块
    score = $3
    species_count = 0
}
# 捕获s开头的行，提取染色体号、起始位置、比对长度，并计数物种数
/^s/ {
    if (species_count == 0) {
        chr = $2
        start = $3
        len = $4  # 使用len代替length
    }
    species_count++
}
# 输出最后一个比对块
END {
    if (species_count > 0) {
        print chr, start, score, len, species_count
    }
}
' chr3.maf | sed 's/ref.chr3/chr3/g' | paste - xx > chr3.count.txt

#######算一下比对的质量_chr4
grep "^s" chr1.maf | awk '{print $2, $3, $4, $5}' | grep "ref" | less -S
grep "score=" chr1.maf | cut -d"=" -f2 |  cut -d"." -f1 > xx
#!/bin/bash
# 输出表头
echo -e "Chromosome\tStart\tScore\tLength\tSpecies_Count"
# 使用awk解析MAF文件
awk '
BEGIN {
    FS = "[ ]+"
    OFS = "\t"
}
# 捕获a score行，并存储比对块的评分
/^a score/ {
    # 如果已经有了一个比对块，则输出它
    if (species_count > 0) {
        print chr, start, score, len, species_count
    }
    # 处理新的比对块
    score = $3
    species_count = 0
}
# 捕获s开头的行，提取染色体号、起始位置、比对长度，并计数物种数
/^s/ {
    if (species_count == 0) {
        chr = $2
        start = $3
        len = $4  # 使用len代替length
    }
    species_count++
}
# 输出最后一个比对块
END {
    if (species_count > 0) {
        print chr, start, score, len, species_count
    }
}
' chr1.maf | sed 's/ref.chr1/chr1/g' | paste - xx > chr1.count.txt




