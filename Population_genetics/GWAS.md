# 1.合并vcf文件（VCF -> VCF）
for file in /public1/home/yinhang/licui/Zhou_data/outfiles/indel_filtered/chr*_indel_filtered.vcf.gz; do
  bcftools index $file
done

bcftools concat -o /public1/home/yinhang/licui/Zhou_data/outfiles/merge_file/merged.vcf.gz -O z /public1/home/yinhang/licui/Zhou_data/outfiles/indel_filtered/chr{1..12}_indel_filtered.vcf.gz

# 2. 拆分为snp和indel的vcf文件（VCF -> vcf）
 # snp
bcftools view -v snps /public1/home/yinhang/licui/Zhou_data/outfiles/10_merge_file/merged.vcf.gz -Oz -o /public1/home/yinhang/licui/Zhou_data/outfiles/10_merge_file/merged_snps.vcf.gz
 # indel
bcftools view -v indels /public1/home/yinhang/licui/Zhou_data/outfiles/10_merge_file/merged.vcf.gz -Oz -o /public1/home/yinhang/licui/Zhou_data/outfiles/10_merge_file/merged_indels.vcf.gz
# 3.plink将vcf转为PLINK格式（VCF -> BIM、BED、FAM）
plink --vcf /public1/home/yinhang/licui/Zhou_data/outfiles/10_merge_file/merged_snps.vcf.gz --make-bed --out /public1/home/yinhang/licui/Zhou_data/outfiles/11_GCTA_data/01_gcta_input_file/SNP/rice_SNP

plink --vcf /public1/home/yinhang/licui/Zhou_data/outfiles/10_merge_file/merged_indels.vcf.gz --make-bed --out /public1/home/yinhang/licui/Zhou_data/outfiles/11_GCTA_data/01_gcta_input_file/INDEL/rice_INDEL

# 4.利用plink做PCA分析(BIM、BED、FAM - > eigenval、eigenvec)
plink --bfile /public1/home/yinhang/licui/Zhou_data/outfiles/11_GCTA_data/01_gcta_input_file/SNP/rice_SNP --pca --out /public1/home/yinhang/licui/Zhou_data/outfiles/11_GCTA_data/01_gcta_input_file/PCA/rice_SNP_pca

plink --bfile /public1/home/yinhang/licui/Zhou_data/outfiles/11_GCTA_data/01_gcta_input_file/INDEL/rice_INDEL --pca --out /public1/home/yinhang/licui/Zhou_data/outfiles/11_GCTA_data/01_gcta_input_file/PCA/rice_INDEL_pca
![](PCA.r) 
![image.png|350](https://image-hosting-website.oss-cn-beijing.aliyuncs.com/img/202310241333522.png)

# 5.GCTA分析（注意-out会生成带有.mlma，所以输入文件不要带有.mlma）
 # 输入文件
 
 # 1.1 计算GRM为GCTA可输入文件(BIM、BED、FAM - > .grm.bin、.grm.id、.grm.N.bin)
 /public1/home/yinhang/software/apps/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile /public1/home/yinhang/licui/Zhou_data/outfiles/11_GCTA_data/01_gcta_input_file/SNP/rice_SNP --make-grm --out /public1/home/yinhang/licui/Zhou_data/outfiles/11_GCTA_data/02_GRM_data/rice_SNP --thread-num 30
 
  /public1/home/yinhang/software/apps/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile /public1/home/yinhang/licui/Zhou_data/outfiles/11_GCTA_data/01_gcta_input_file/INDEL/rice_INDEL --make-grm --out /public1/home/yinhang/licui/Zhou_data/outfiles/11_GCTA_data/02_GRM_data/rice_INDEL --thread-num 30
  
  # 1.2 确定PCA数量后，整理PCA结果文件为GCTA可输入文件（eigenvec - > result.txt）
![](PCA.r)
```R
#生成特定格式PCA结果
eigenvec <- read.table("D:/Desktop/rice_SNP_pca.eigenvec", header=FALSE)
# 选择前3个主成分和样本ID
# 样本ID在第1列，主成分在第2, 3, 4列
pca_result <- eigenvec[, c(1, 2, 3, 4, 5)]  
# 添加列名
colnames(pca_result) <- c("Family_ID","Iyellowividual_ID","PC1", "PC2", "PC3")
write.table(pca_result, "D:/Desktop/rice_SNP_pca_result_3.txt", sep="\t", row.names=FALSE, quote=FALSE)
```
# 1.3 整理表型文件，为GCTA可输入文件（phenotye.txt - > result.txt）

```shell
  # 循环处理每个表型文件
for file in "${phenotype_dir}"*.txt; do
    # 提取文件名（不包含路径和扩展名）
    filename=$(basename -- "$file")
    filename_noext="${filename%.*}"
    # 生成Phenotype Plink格式输出文件路径
    plink_output_file="${output_dir}03_phenotype_plink/${filename_noext}.txt"
    # 执行表型文件转换为Plink格式命令
    awk 'NR>1 {print $1, $1, $2}' "$file" > "$plink_output_file"
#awk 'NR>1 {print $1, $1, $2}' "/public1/home/yinhang/projects/licui/Zhou_data/04_outfiles/11_GCTA_data/03_phenotype_plink/XQ.txt" > "/public1/home/yinhang/projects/licui/Zhou_data/04_outfiles/11_GCTA_data/03_phenotype_plink/XQ1.txt"
```
# 输出文件 
  # 1.脚本Heritability_GWAS（PLINK格式文件、GRM、表型文件、PCA - > .mlma）
```
    for grm_type in "_SNP" "_INDEL"; do
        # 生成遗传率结果输出文件路径
        heritability_output_file="${output_dir}04_Heritability/${filename_noext}${grm_type}_heritability.log"
        # 生成GWAS结果输出文件路径
        gwas_output_file="${output_dir}05_GWAS_results/${filename_noext}${grm_type}_gwas.log"
        # 执行遗传率分析命令
        /public1/home/yinhang/software/apps/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm "${output_dir}02_GRM_data/rice${grm_type}" --pheno "$plink_output_file" --reml --out "$heritability_output_file" --thread-num 30
        # 执行GWAS分析命令for SNP
        /public1/home/yinhang/software/apps/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile "${output_dir}01_gcta_input_file/SNP/rice_SNP" --grm "${output_dir}02_GRM_data/rice_SNP${grm_type}" --pheno "$plink_output_file" --qcovar "/public1/home/yinhang/licui/Zhou_data/outfiles/11_GCTA_data/01_gcta_input_file/PCA/pca_result.txt" --mlma --out "$gwas_output_file" --thread-num 30