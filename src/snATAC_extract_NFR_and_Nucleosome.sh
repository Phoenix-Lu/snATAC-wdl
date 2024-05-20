#! /bin/bash
bed_file=$1

base_name=$(basename "${bed_file}")
base_name=${base_name%.*}

absolute() {
    if (( $1 < 0 )); then
        echo $(( -$1 ))
    else
        echo $1
    fi
}

# 打开输出文件
touch "${base_name}_Nucleosome.bed" "${base_name}_NFR.bed"

# 循环读取输入文件，并将记录写入对应的输出文件
# while IFS=$'\t' read -r col1 col2 col3 col4; do
#     # 计算第二列和第三列的绝对值
#     diff=$(absolute $((col3 - col2)))
#     if (( diff > 100 )); then
#         echo "$col1 $col2 $col3 $col4" >> "${base_name}_Nucleosome.bed"
#     else
#         echo "$col1 $col2 $col3 $col4" >> "${base_name}_NFR.bed"
#     fi
# done < "${bed_file}"

awk 'BEGIN{OFS="\t"} {
    diff = ($3 - $2) < 0 ? -($3 - $2) : ($3 - $2);
    if (diff > 100) {
        print $0 >> "'${base_name}_Nucleosome.bed'";
    } else {
        print $0 >> "'${base_name}_NFR.bed'";
    }
}' "${bed_file}"
# 关闭输出文件
