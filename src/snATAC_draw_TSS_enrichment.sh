#!/bin/bash

bw_file=$1
gtf_filt=$2


#TODO 需要测试，暂时缺文件，需要完善，有可能出黑色是这个的问题
# gtf转bed
sed 's/"/\t/g' "${gtf_filt}" | \
awk 'BEGIN{OFS=FS="\t"}{if($3=="gene") {if($7=="+") {start=$4-1000; end=$4+500;} else {if($7=="-") start=$5-500; end=$5+1000; } if(start<0) start=0; print $1,start,end,$14,$10,$7;}}' > computeMtrix.generef
computeMatrix reference-point  --referencePoint TSS  -p 15 -b 3000 -a 3000 -R computeMtrix.generef -S "${bw_file}" --skipZeros -o matrix_test_TSS.gz --outFileSortedRegions regions_test_genes.bed

plotHeatmap -m matrix_test_TSS.gz -out test_Heatmap.png

# reference-point # 选择模式
# -p 15 线程
# --referencePoint TSS  # 选择参考点: TES, center
# -b 10000 -a 10000  # 感兴趣的区域，-b上游，-a下游
# -R  基因注释信息
# -S  提供的 bigwig 文件
# --skipZeros 
# -out ./test.TSS.gz  输出为文件用于plotHeatmap, plotProfile
#--outFileSortedRegions  ./test.TSS.bed  输出的文件名