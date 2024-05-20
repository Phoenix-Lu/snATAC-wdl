#!/bin/bash



### 基因组注释
blacklist_ref="/work/xulab/xulab-seq/reference/annotation/encode/blacklist/hg38-blacklist.v2.bed"
nucleosome_ref="/work/xulab/xulab-seq/reference/annotation/paper/nucleosome_position/merge.nucleosome_pos_130_160.bed"

PLS_ref="/work/xulab/xulab-seq/reference/annotation/encode/SCREEN/hg38/ALL/GRCh38-PLS.bed"
pELS_ref="/work/xulab/xulab-seq/reference/annotation/encode/SCREEN/hg38/ALL/GRCh38-pELS.bed"
dELS_ref="/work/xulab/xulab-seq/reference/annotation/encode/SCREEN/hg38/ALL/GRCh38-dELS.bed"
ctcf_ref="/work/xulab/xulab-seq/reference/annotation/encode/SCREEN/hg38/ALL/GRCh38-CTCF.bed"
genebody_ref="/work/xulab/xulab-seq/reference/annotation/gencode/hg38/genebody_excluded_tss.bed"

ls ${bed_dir}/*.bed.gz | grep -iE 'NFR|Nucleosome' | grep -i -v pi | while read id;
do

file=$(basename $id)
sample=${file/.bed.gz}

# 定义输出文件
output_csv="${bed_dir}/${sample}_counts.csv"

# 初始化CSV文件并写入表头
echo "Sample,TotalReads,Blacklist,Promoter,pEnhancer,dEnhancer,CTCF,Genebody,Nucleosome" > "$output_csv"

# 初始化包含样本名称的行
line="${sample}"
# 计算总的reads数量
totalReads=$(zcat "${id}" | wc -l)
line="${line},${totalReads}"

for ref in blacklist_ref PLS_ref pELS_ref dELS_ref ctcf_ref genebody_ref nucleosome_ref;
do

region=${ref/_ref/}

ref_file=${!ref}

# 计算在当前区域内的reads数量
regionReads=$(bedtools sort -i "${ref_file}" | bedtools merge -i stdin | bedtools intersect -u -nonamecheck -a "${id}" -b stdin | wc -l)

# 将结果添加到行中
line="${line},${regionReads}"

done

# 写出CSV文件
echo "${line}" >> "$output_csv"
done