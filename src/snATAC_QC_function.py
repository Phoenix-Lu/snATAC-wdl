import subprocess
import os
import argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os
import pandas as pd
from jinja2 import Environment, FileSystemLoader
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import sys
import jinja2
import base64
import re


def fastqc_report():
    pass

def classify_fragment(length):
    if length < 100:
        return 'NFR'
    elif 100 <= length < 300:
        return 'Mono_nucle'
    elif 300 <= length < 500:
        return 'Dimer_nucle'
    elif 500 <= length < 700:
        return 'Tri_nucle'
    elif 700 <= length < 900:
        return 'Tetra_nucle'
    else:
        return 'Multi_nucle'

# TESTPASSED
def get_chrM_ratio(input_bam, pattern = 'chrM'):
    all_reads_num = subprocess.run(f"samtools view {input_bam} | wc -l ",
                                   shell=True, capture_output=True,encoding='utf-8').stdout
    print(all_reads_num)
    chrm_reads_num = subprocess.run(f"samtools view {input_bam} | grep {pattern} | wc -l",
                                    shell=True, capture_output= True,encoding='utf-8').stdout
    chrm_ratio = int(chrm_reads_num)/int(all_reads_num)
    print(chrm_reads_num)
    file = './QCs/' + os.path.basename(input_bam).split('.')[0] + '.chrmratio'
    with open(file, 'a') as file:
        print(f"chrm_ratio for {input_bam} is {chrm_ratio}")
        print(chrm_ratio, file=file)
    return chrm_ratio

# TESTPASSED
def get_bam_qc_report(input_bam, filter_para=None):
    qc_report = os.path.basename(input_bam).split('.')[0] + ".filtersummary"
    # 打印stat报告
    with open(qc_report, "a") as f:
        print(
            f"{'*' * 20} {input_bam} full mapping QC Report by samtools stat {'*' * 20}",
            file=f,
        )
        print(filter_para, file=f)
    command = f"samtools stat {input_bam} |grep ^SN | cut -f 2- >> {qc_report}"
    subprocess.run(command, shell=True)
    # 打印flagstat报告
    with open(qc_report, "a") as f:
        print("\n\n", file=f)
        print(
            f"{'*' * 20} {input_bam} full mapping QC Report by samtools flagstat {'*' * 20}",
            file=f,
        )
    command = f"samtools flagstat {input_bam} >> {qc_report}"
    subprocess.run(command, shell=True)
    return qc_report

# TESTPASSED
def get_insert_length_distribution(input_bam):
    import os
    import subprocess
    base_name = os.path.basename(input_bam).split(".")[0]
    insert_len_file = './QCs/'+ base_name + ".insertlen"
    # 双端变单端，获得的是每次测序的insert length
    command = f"""samtools view {input_bam} | awk -F'\\t' 'function abs(x){{return ((x < 0.0) ? -x : x)}} {{print $1"\\t"abs($9)}}' | sort | uniq > {insert_len_file}"""
    subprocess.run(command, shell=True)
    len_distribute = './QCs/'+base_name + ".lenDistribute"
    command2 = f"""awk '{{print $2}}' {insert_len_file} | sort | uniq -c | awk '{{print $2 "\t" $1}}' > {len_distribute}"""
    subprocess.run(command2, shell=True)
    return insert_len_file, len_distribute


# unused
def draw_insert_length_distribution_plot(insert_len_file):
    import pandas as pd
    import matplotlib.pyplot as plt
    import os
    base_name = os.path.basename(insert_len_file).split(".")[0]
    data = pd.read_csv(insert_len_file, sep="\t", header=None)
    data.columns = ["read_name", "insert_length"]
    data = data[data["insert_length"] < 2000]
    sizes = list(data["insert_length"].value_counts())
    plt.hist(sizes, bins=500)
    plt.xlabel("Insert Length")
    plt.ylabel("Read Counts")
    plt.title(f"Insert Length Distribution for {base_name}")
    plt.savefig(base_name + ".png")
    return None


# TESTPASSED
def get_align_qc_report(input_bam, filter_para=None):
    basename = (os.path.basename(input_bam)).split('.')[0]
    align_qc_report_stat = './QCs/' + basename + ".stat_summary"
    # 打印stat 报告
    with open(align_qc_report_stat, "a") as f:
        print(
            f"{'*' * 20} {input_bam} full mapping QC Report by samtools stat {'*' * 20}",
            file=f,
        )
        print(filter_para, file=f)
    command = f"samtools stat {input_bam} |grep ^SN | cut -f 2- >> {align_qc_report_stat}"
    subprocess.run(command, shell=True)
    align_qc_report_flagstat = './QCs/' + basename + ".flagstate_summary"
    # 打印flagstat报告
    with open(align_qc_report_flagstat, "a") as f:
        print("\n\n", file=f)
        print(
            f"{'*' * 20} {input_bam} full mapping QC Report by samtools flagstat {'*' * 20}",
            file=f,
        )
    command = f"samtools flagstat {input_bam} >> {align_qc_report_flagstat}"
    subprocess.run(command, shell=True)
    return align_qc_report_stat, align_qc_report_flagstat


def get_align_qc_report_summary(full_qc_report, prefix, sample_name):
    align_qc_report = sample_name + prefix + ".alignsummary"




# TESTPASSED
def calculate_region_enrichment_score(reads_bed, region_bed):
    import os
    import subprocess
    command = "bedtools sort -i {} | bedtools merge -i stdin | bedtools intersect -u -nonamecheck -a {} -b stdin | wc -l".format(region_bed, reads_bed)
    intersect_read_count = int(subprocess.run(command, shell=True, capture_output=True, encoding="utf-8").stdout)
    total_read_count = int(subprocess.run(f"wc -l {reads_bed}", shell=True, capture_output=True, encoding="utf-8").stdout.split()[0])
    # print(intersect_read_count)
    # print(total_read_count)
    fract_reads = float(intersect_read_count) / total_read_count
    return fract_reads

# TESTPASSED
def get_gf_annotation_csv(reads_beds, region_ref, sample_name):
    gf_annotation_csv = []
    for reads_bed in reads_beds:
        basename = (os.path.basename(reads_bed)).split('.')[0]
        header = []
        ratio = []
        with open(reads_bed,'r') as file:
            reads_num = len(file.readlines())
        for ref_name, ref_file in region_ref.items():
            header.append(ref_name.split('_')[0])
            ratio.append(calculate_region_enrichment_score(reads_bed,ref_file))
        # 记录的是数字，而不是ratio
        csv = pd.DataFrame([[basename + '_counts', reads_num] + [int(fract * reads_num) for fract in ratio], [basename + '_portion', '1.0'] + ratio], columns = ['Sample', 'TotalReads'] + header)
        gf_annotation_csv.append(csv)
    gf_annotation_csv = pd.concat(gf_annotation_csv, axis = 0, ignore_index = True)
    gf_annotation_csv = gf_annotation_csv.T
    gf_annotation_csv = gf_annotation_csv.reset_index(drop=False)
    gf_annotation_csv.columns = range(gf_annotation_csv.shape[1])

    # 重新将第一行设置为列名
    new_header = gf_annotation_csv.iloc[0]  # 抓取第一行作为header
    gf_annotation_csv = gf_annotation_csv[1:]  # 取数据少掉header的部分
    gf_annotation_csv.columns = new_header  # 设置新的header
    gf_annotation_csv.reset_index(drop=True, inplace=True)  # 重置索引并删除旧的索引列
    # gf_annotation_csv.drop('Sample')
    gf_annotation_csv.to_csv(f'./QCs/{sample_name}_annotation_counts.csv', index=True)
    return gf_annotation_csv

    
def get_gf_frag_len(fragment_counts_path, sample_name):
    ### 片段长度分布
    fragLen = pd.DataFrame()
    file_name = fragment_counts_path
    sampleInfo = sample_name
    # 读取文件
    temp_df = pd.read_table(file_name, header=None, names=['V1', 'V2'])
            
    # 转换数据类型并计算Weight
    temp_df['fragLen'] = temp_df['V1'].astype(float)
    temp_df['fragCount'] = temp_df['V2'].astype(float)
    temp_df['Weight'] = temp_df['V2'] / temp_df['V2'].sum()
    temp_df['sampleInfo'] = sampleInfo
            
    # 合并到最终的DataFrame
    fragLen = pd.concat([fragLen, temp_df], ignore_index=True)

    fragLen_log = fragLen.copy()
    fragLen_log['fragCount'] = np.log10(fragLen_log['fragCount'])

    # 应用这个函数到DataFrame上创建一个新列
    
    fragLen['FragmentType'] = fragLen['fragLen'].apply(classify_fragment)

    # 对每个 sampleInfo 和 FragmentType 组合的 fragCount 进行求和
    fragment_counts_sum = fragLen.groupby(['sampleInfo', 'FragmentType'])['fragCount'].sum().reset_index(name='TotalFragCount')

    # 计算每个样本的总 fragCount
    total_frag_count_per_sample = fragment_counts_sum.groupby('sampleInfo')['TotalFragCount'].sum().reset_index(name='TotalCountPerSample')

    # 合并总数回到原始计数 DataFrame，为了计算比例
    fragment_counts_sum = fragment_counts_sum.merge(total_frag_count_per_sample, on='sampleInfo')

    # 计算比例
    fragment_counts_sum['Proportion'] = fragment_counts_sum['TotalFragCount'] / fragment_counts_sum['TotalCountPerSample']

    # 使用 pivot_table 创建数量和比例的表格
    # 对于数量
    pivot_table_counts = fragment_counts_sum.pivot_table(index='sampleInfo', columns='FragmentType', values='TotalFragCount', fill_value=0)

    # 对于比例
    pivot_table_proportions = fragment_counts_sum.pivot_table(index='sampleInfo', columns='FragmentType', values='Proportion', fill_value=0)

    pivot_table_proportions = pivot_table_proportions.map(lambda x: round(x, 2))

    # 定义列的期望顺序
    desired_order = ['NFR', 'Mono_nucle', 'Dimer_nucle', 'Tri_nucle', 'Tetra_nucle', 'Multi_nucle']

    # 确保 pivot_table_counts 和 pivot_table_proportions 包含所有期望的列，不存在的列设置为0
    for table in [pivot_table_counts, pivot_table_proportions]:
        for col in desired_order:
            if col not in table.columns:
                table[col] = 0  # 为不存在的列添加0

    # 现在所有期望的列都存在了，可以安全地按照 desired_order 进行重排序
    pivot_table_counts = pivot_table_counts[desired_order]
    pivot_table_proportions = pivot_table_proportions[desired_order]

    # 重置索引以便合并
    pivot_table_counts.reset_index(inplace=True)
    pivot_table_proportions.reset_index(inplace=True)

    # 上下合并
    merged_table_fragLen = pd.concat([pivot_table_counts, pivot_table_proportions], ignore_index=True)

    merged_table_fragLen.sort_values(by=['sampleInfo'], inplace=True)


    data = fragLen[fragLen["sampleInfo"] == sampleInfo ]
    data_log = fragLen_log[fragLen_log["sampleInfo"] == sampleInfo ]
    # 设置绘图风格
    sns.set_theme(style="white")
    
    # 创建图形和轴对象
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # 绘制线图
    sns.lineplot(x='fragLen', y='fragCount', hue='sampleInfo', data=data, ax=ax, palette="magma")
    
    # 设置x轴和y轴的标签
    ax.set_xlabel("Fragment Length", fontsize=20)
    ax.set_ylabel("Count", fontsize=20)
    
    # 设置标题
    ax.set_title("Fragment Length Distribution", fontsize=24)
    
    # 设置x轴的范围
    ax.set_xlim(0, 800)
    
    ax.get_legend().remove()

    # 在大图中固定位置添加虚线
    for x in [100, 300, 500, 700, 900]:
        ax.axvline(x=x, color='grey', linestyle='--', linewidth=1)

    sns.despine()
    
    ### 第二张图：在第一张图的右上角添加第二张小图
    axins = inset_axes(ax, width="55%", height="55%", loc='upper right',
                    bbox_to_anchor=(0.29, 0.29, 0.7, 0.7),
                    bbox_transform=ax.transAxes)
    
    # 在小图中绘制对数变换后的数据
    sns.lineplot(x='fragLen', y='fragCount', hue='sampleInfo', data=data_log, ax=axins, palette="magma")
    axins.set_title("Log Scale", fontsize=10)
    axins.set_xlabel('')
    axins.set_ylabel('')
    # 设置小图的坐标轴范围，可能需要根据对数变换后的数据进行调整
    axins.set_xlim(0, 800)
    
    # 设置y轴的范围以匹配对数变换的值
    # axins.set_ylim(data_log['fragCount'].min(), fragLen_log['fragCount'].max() + 1)
    
    sns.despine(ax = axins)
    
    # 保存图像到文件，指定分辨率和格式
    out_path_png = './QCs/' + sampleInfo+ '_fragment_length_plot.png'
    plt.savefig(out_path_png, dpi=300, format='png')

    df = merged_table_fragLen[merged_table_fragLen['sampleInfo'] == sampleInfo].iloc[:,1:]

    df.insert(0, 'Type', ['Counts', 'Proportions'])

    df.iloc[0:1,2:] = df.iloc[0:1,2:].astype(int).astype(str)
    out_path_csv = './QCs/' + sampleInfo + '_fragment_length.csv'

    df.to_csv(out_path_csv)
    return out_path_png, df


# TESTPASSED
def get_gf_mapAndRmdup_csv(align_result, states_result, chrom_ratio, accept_states_result, sample, picard_metric):
    with open(align_result, 'r') as f:
        align_result_list = f.readlines()
    with open(states_result, 'r') as f:
        states_result_list = f.readlines()
    with open(accept_states_result, 'r') as f:
        accept_states_result_list = f.readlines()
    start_index = [i for i, line in enumerate(align_result_list) if 'reads; of these:' in line][0]
    total_reads = int(align_result_list[start_index].split()[0])
    # 初始化计数器
    aligned_zero = 0
    aligned_one = 0
    aligned_more = 0
    # 解析每个段落
    aligned_zero += int(align_result_list[start_index + 2].split()[0])
    aligned_one += int(align_result_list[start_index + 2 +1].split()[0])
    aligned_more += int(align_result_list[start_index + 2 + 2].split()[0])
    overall_alignment_rate = round((aligned_one + aligned_more) / total_reads * 100, 2)
        # 计算未比对上的reads数
    aligned_zero = total_reads - (aligned_one + aligned_more)
    chrM_mapped_ratio = chrom_ratio
    Exact_mapped_reads = int(int(states_result_list[7].split('\t')[-1]) / 2)
    Exact_mapped_reads_rate = round(int(states_result_list[7].split('\t')[-1]) / 2 / total_reads * 100, 2)
    Exact_rmdup_mapped_reads = int(int(accept_states_result_list[7].split('\t')[-1]) / 2)
    Exact_rmdup_mapped_reads_rate = round(int(accept_states_result_list[7].split('\t')[-1]) / 2 / total_reads * 100, 2)
    align_dataframe = pd.DataFrame({
    'Sample': [sample],
    'TotalReads': [total_reads],
    'AlignmentReads': [aligned_one + aligned_more],
    'AlignmentRate': [overall_alignment_rate],
    'ChrM_mapped_ratio':[chrM_mapped_ratio],
    'ExactMappedReads': [Exact_mapped_reads],
    'ExactMappedRate': [Exact_mapped_reads_rate],
    'ExactMappedReads_rmdup': [Exact_rmdup_mapped_reads],
    'ExactMappedRate_rmdup': [Exact_rmdup_mapped_reads_rate]
})
    ### PCR duplitcates
        # 初始化一个空的DataFrame来存储最终结果
    dupResult = pd.DataFrame()
    file_name = picard_metric
    # 读取文件
    dupRes = pd.read_csv(file_name, sep='\t', comment='#', header=0)
    dupRes.iloc[2:, 2] = pd.to_numeric(dupRes.iloc[2:, 2], errors='coerce')
    single = dupRes.iloc[2,2]
    n_duplicate_set = dupRes.iloc[2:, 2].sum()
    Saturation = round((1 - single / n_duplicate_set) * 100, 2)
    # 创建一个临时DataFrame存储计算结果
    temp_df = pd.DataFrame({
    'Sample': [sample],
    'MappedFragNum': [int(dupRes['READ_PAIRS_EXAMINED'].iloc[0])],
    'DuplicationRate': [round(float(dupRes['PERCENT_DUPLICATION'].iloc[0]) * 100, 2)],
    'EstimatedLibrarySize': [int(dupRes['ESTIMATED_LIBRARY_SIZE'].iloc[0])],
    'SequencingSaturation': [Saturation],
    'UniqueFragNum': [n_duplicate_set]
    })
    # 合并到最终的DataFrame
    mapAndRmdup = pd.concat([align_dataframe, temp_df],axis = 1)
    mapAndRmdup = mapAndRmdup.T
    mapAndRmdup = mapAndRmdup.reset_index(drop=False)
    mapAndRmdup.columns = range(mapAndRmdup.shape[1])

    # 重新将第一行设置为列名
    new_header = mapAndRmdup.iloc[0]  # 抓取第一行作为header
    mapAndRmdup = mapAndRmdup[1:]  # 取数据少掉header的部分
    mapAndRmdup.columns = new_header  # 设置新的header
    mapAndRmdup.reset_index(drop=True, inplace=True)  # 重置索引并删除旧的索引列
    return mapAndRmdup

def get_filter_qc_report(input_bam, prefix, sample_name, filter_para=None):
    filter_qc_report = sample_name + prefix + ".filtersummary"
    # 打印stat报告
    with open(filter_qc_report, "a") as f:
        print(
            f"{'*' * 20} {input_bam} full mapping QC Report by samtools stat {'*' * 20}",
            file=f,
        )
        print(filter_para, file=f)
    command = f"samtools stat {input_bam} |grep ^SN | cut -f 2- >> {filter_qc_report}"
    subprocess.run(command, shell=True)

    # 打印flagstat报告
    with open(filter_qc_report, "a") as f:
        print("\n\n", file=f)
        print(
            f"{'*' * 20} {input_bam} full mapping QC Report by samtools flagstat {'*' * 20}",
            file=f,
        )
    command = f"samtools flagstat {input_bam} >> {filter_qc_report}"
    subprocess.run(command, shell=True)
    return filter_qc_report
    

def gf_seq_saturation(metric, sample):
    import os
    import sys
    import random
    import argparse
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    from collections import Counter

    def generate_elements(bins, counts):
        n = 0
        elements = []
        rev_bins, rev_counts = bins[::-1], counts[::-1]
        for idx, bin in enumerate(rev_bins):
            count = rev_counts[idx]
            for _ in range(count):
                for _ in range(bin):
                    elements.append(n)
                n+=1
        return elements


    def calculate_bin_counts(elements):
        if not elements:
            return [], []
        element_counts = Counter(elements)

        value_counts = {}
        for key, value in element_counts.items():
            if value in value_counts:
                value_counts[value] += 1
            else:
                value_counts[value] = 1
        sorted_keys = sorted(value_counts.keys())
        return sorted_keys, [value_counts[key] for key in sorted_keys]


    def sample_bins_and_counts(elements, sample_ratio):
        sample_size = int(len(elements) * sample_ratio)
        sample_elements = random.sample(elements, sample_size)
        new_bins_transformed, new_counts_transformed = calculate_bin_counts(sorted(sample_elements))
        
        return new_bins_transformed, new_counts_transformed


    def compute_seq_saturation(bins, counts):
        single = counts[0]
        n_duplicate_set = sum(counts)
        total, dup = 0, 0
        for i in range(0, len(bins)):
            total += bins[i]*counts[i]
            dup += (bins[i]-1)*counts[i]
        PERCENT_DUPLICATION = (dup)*100/(total)
        seq_saturation = (1-(single/n_duplicate_set))*100
        return round(seq_saturation, 2), single, n_duplicate_set, round(PERCENT_DUPLICATION, 2)


    def write_tuple_to_file(file_handle, rate, tuple_data):
        tuple_str = '\t'.join(map(str, tuple_data))
        file_handle.write(str(rate) + "\t" + tuple_str + '\n')

    metric = os.path.abspath(metric)
    out_plot_par = os.path.abspath(os.path.join('./QCs', sample+ 'seq_saturation.bar'))

    df = pd.read_csv(metric, sep="\t", skiprows=11, header=None)
    bins = [int(num) for num in list(df[0])]
    counts = [int(num) for num in list(df[2])]

    elements = generate_elements(bins, counts)
    print(f"finish generate_elements")

    f = open(out_plot_par, "w")
    for num in [0.0001, 0.001, 0.01, 0.1]:
        for i in range(1, 10):
            sample_ratio = round(num*i, 4)
            new_bins_transformed, new_counts_transformed = sample_bins_and_counts(elements, sample_ratio)
            write_tuple_to_file(f, sample_ratio, compute_seq_saturation(new_bins_transformed, new_counts_transformed))
            print(f"finish {sample_ratio}")
    write_tuple_to_file(f, 1, compute_seq_saturation(bins, counts))
    f.close()
    return out_plot_par
def gf_plot_seq_saturation(bar_file_path, sample):
    import os
    import sys
    import argparse
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd


    file_path = os.path.abspath(bar_file_path)
    file_dir = './QCs/'

    df = pd.read_csv(file_path, sep="\s+", header=None)

    df.columns = ['name', 'sequencing_saturation', 'n_single_dupset', 'n_multi_dupset', 'dup_rate']

    print(df.head())


    def plot_sequencing_saturation(df, x_end, out):
        # 绘图
        plt.figure(figsize=(10, 6))

        plt.plot(df['name'], df['sequencing_saturation'], label='Sequencing Saturation', color='blue', marker='o', linestyle='-', linewidth=2, markersize=1)
        plt.plot(df['name'], df['dup_rate'], label='Duplication Rate', color='red', marker='o', linestyle='-', linewidth=2, markersize=1)

        # 设置坐标轴范围和标签
        plt.xlim(0, x_end)
        # plt.xlim(0, 0.08)
        plt.ylim(0, 100)
        plt.xlabel('Downsample Ratio', fontsize=18)  # 横轴标题
        plt.ylabel('Ratio', fontsize=18)  # 纵轴标题

        # 添加图例和网格
        plt.legend()
        plt.grid(True)

        # 保存为PDF
        plt.savefig(out, format='png')
        # plt.savefig('sequencing_saturation.pdf', format='pdf')
    plot_sequencing_saturation(df, 1, os.path.join(file_dir, f"{sample}_sequencing_saturation_1.png"))
    return os.path.join(file_dir, f"{sample}_sequencing_saturation_1.png")

# 提取分析结果中的各项指标并整理成分析报告的形式
def gf_html_gen(sample_name, qc_result_csv, align_result_csv, saturation_image, frag_len_image, fraglen_csv, annotation_csv):
    import os
    import pandas as pd
    from jinja2 import Environment, FileSystemLoader
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    import sys
    import jinja2
    import base64
    template_dir = '/work/xulab/xulab-seq/gf/templates'
    template_file = 'template.html'

    # 初始化Jinja2环境
    env = Environment(loader=FileSystemLoader(template_dir))
    template = env.get_template(template_file)

    # 设置标题为工作目录的最后一个部分
    title = sample_name

    def image_to_base64(image_path):
        with open(image_path, 'rb') as image_file:
            encoded_string = base64.b64encode(image_file.read()).decode('utf-8')
        return encoded_string

    try:
        result_dir ="./QCs"
        os.makedirs(result_dir, exist_ok = True)

        ### 第一部分：质控结果
        # 将转置后的DataFrame转换为可直接在HTML模板中遍历的记录列表
        qc_result = qc_result_csv.to_dict(orient='records')
        
        ### 第二部分：比对结果 
        
        alignResult_sample = align_result_csv.to_dict(orient='records')

        # 测序饱和度图
        saturation_image_path = saturation_image

        # 检查图片是否存在并转换为Base64编码
        if os.path.isfile(saturation_image_path):
            saturation_image_base64 = image_to_base64(saturation_image_path)
        else:
            saturation_image_base64 = None

        ### 第三部分：片段长度分布
        fragLen_image_path = frag_len_image

        # 检查图片是否存在并转换为Base64编码
        if os.path.isfile(fragLen_image_path):
            fragLen_image_base64 = image_to_base64(fragLen_image_path)
        else:
            fragLen_image_base64 = None

        # 片段长度分布表
        fragLen_sample = fraglen_csv.to_dict(orient='records')
        ### 第四部分：注释基因组区域丰度
        genomic_region_annotation = annotation_csv.to_dict(orient='records')
        sample = [{'id': sample_name, 
                        'name': sample_name, 
                        'qc_data': qc_result, 
                        'align_result':alignResult_sample, 
                        'saturation_image_base64':saturation_image_base64,
                        'fragLen_image_base64':fragLen_image_base64,
                        'fragLen_result':fragLen_sample,
                        'genomic_region_annotation':genomic_region_annotation
                    }]
            
        # 渲染HTML
        rendered_html = template.render(title=title, samples=sample)

        out_path = result_dir + '/' + sample_name + '数据分析报告.html'
        # 输出结果到HTML文件
        with open(out_path, 'w') as f:
            f.write(rendered_html)

    except:
        import sys
        print (sys.exc_info()[0])
        import traceback
        print (traceback.format_exc())
    finally:
        print("\n")
        print("Finished !!!")

def parse_arguments():
    parser = argparse.ArgumentParser(prog="Phoenix macs2 call peak.")
    parser.add_argument('--no_filter_bam_file')
    parser.add_argument('--accept_bam_file')
    parser.add_argument('--align_summary')
    parser.add_argument('--picard_metric')
    parser.add_argument('--sample_name', required=True)
    parser.add_argument('--full_bed_file')
    parser.add_argument('--part_bed_files', nargs='+')
    args = parser.parse_args()
    return args


if __name__=="__main__":
    args = parse_arguments()
    try:
        os.mkdir('./QCs')
    except:
        pass
    REGION_REF = {
    'blacklist_ref':"/work/xulab/xulab-seq/reference/annotation/encode/blacklist/hg38-blacklist.v2.bed",
    'nucleosome_ref':"/work/xulab/xulab-seq/reference/annotation/paper/nucleosome_position/merge.nucleosome_pos_130_160.bed",
    'PLS_ref':"/work/xulab/xulab-seq/reference/annotation/encode/SCREEN/hg38/ALL/GRCh38-PLS.bed",
    'pELS_ref':"/work/xulab/xulab-seq/reference/annotation/encode/SCREEN/hg38/ALL/GRCh38-pELS.bed",
    'dELS_ref':"/work/xulab/xulab-seq/reference/annotation/encode/SCREEN/hg38/ALL/GRCh38-dELS.bed",
    'ctcf_ref':"/work/xulab/xulab-seq/reference/annotation/encode/SCREEN/hg38/ALL/GRCh38-CTCF.bed",
    'genebody_ref':"/work/xulab/xulab-seq/reference/annotation/gencode/hg38/genebody_excluded_tss.bed"
    }
    align_result= args.align_summary
    bam_file = args.no_filter_bam_file
    accept_bam_file = args.accept_bam_file
    picard_metric = args.picard_metric
    full_bed_file = args.full_bed_file
    parts_bed_file = args.part_bed_files
    sample_name = args.sample_name
    if (bam_file):
        chrom_ratio = get_chrM_ratio(bam_file)
        bam_state = get_align_qc_report(bam_file)[0]
    else:
        print("*******chrom_ratio, bam_state\n cannot be calculated because bam_without_filter is not provided, chrom_ratio is set to -1")
        chrom_ratio = -1

    if (accept_bam_file):
        accept_bam_state = get_align_qc_report(accept_bam_file)[0]
        insert_len_file = get_insert_length_distribution(accept_bam_file)[1]
        frag_len_image, fraglen_csv = get_gf_frag_len(insert_len_file, sample_name)
    else:
        print("*******accept_bam_state,insert_len_file,frag_len_image,fraglen_csv\n cannot be calculated because accept_bam_file is not provided")
    
    if (picard_metric):
        saturation_bar = gf_seq_saturation(picard_metric, sample=sample_name)
        saturation_bar_plot = gf_plot_seq_saturation(saturation_bar, sample = sample_name)
    else:
        print("*******picard_metric not provided, saturation is not calculated.")


    if (accept_bam_file and bam_file and picard_metric and align_result):
        align_result_csv = get_gf_mapAndRmdup_csv(align_result, bam_state, chrom_ratio, accept_bam_state,sample_name , picard_metric)
    else:
        print("*******align_result_csv is not caculated")

    if (parts_bed_file or full_bed_file):
        if parts_bed_file == None:
            parts_bed_file = []
        if full_bed_file == None:
            full_bed_file = []
        else:
            full_bed_file = [full_bed_file]
        annotation_csv = get_gf_annotation_csv(full_bed_file + parts_bed_file,REGION_REF,sample_name)
    else:
        print("*******no bed file provided, no annotation_csv")
    qc_result_csv = pd.DataFrame({'':["QC_report", "See fastp report instead of here"]})
    
    if(align_result and bam_file and accept_bam_file and picard_metric and (full_bed_file or parts_bed_file)):
        gf_html_gen(sample_name, qc_result_csv, align_result_csv, saturation_bar_plot, frag_len_image, fraglen_csv, annotation_csv)
    else:
        print("*********no html generated")