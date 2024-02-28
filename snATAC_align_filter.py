# Author: Phoenix_Lu

import sys
import os
import argparse
import subprocess


def parse_arguments():
    parser = argparse.ArgumentParser(prog="Phoenix Bam filter.")
    parser.add_argument("--bam_file", required=True)
    parser.add_argument("--mpQ", default="30", type=str)
    parser.add_argument("--rmunpro_pair", default=False, action="store_true")
    parser.add_argument("--rmunmap", default=False, action="store_true")
    parser.add_argument("--rmdup", default=False, action="store_true")
    parser.add_argument("--rmchrM", default=False, action="store_true")
    parser.add_argument("--threads", default=5)
    parser.add_argument("--prefix", default="_filtered", type=str)
    parser.add_argument("--sample_name", required=True)
    args = parser.parse_args()
    return args


def do_bam_filter(
    bam_file, prefix, rmdup, mpQ_threhold, rmUnpro_pair, rmUnmap, rmchrM, sample_name
):

    output_bam = sample_name + prefix + ".bam"
    # para_prefix是为了后续报告提供的
    para_prefix = {
        "mapQ": mpQ_threhold,
        "remove duplicates": True if rmdup else False,
        "remove Unproperly paired": True if rmUnpro_pair else False,
        "remove unmaped": True if rmUnmap else False,
        "remove chrM": True if rmchrM else False,
    }

    # 加入command参数
    command = ["samtools", "view", bam_file, "-h"]

    # 加入参数
    if int(mpQ_threhold) > 0:
        command.extend(["-q", f"{mpQ_threhold}"])

    if rmdup or rmUnmap:
        command.extend(["-F", ("1024," if rmdup else "") + ("4," if rmUnmap else "")])

    if rmUnpro_pair:
        command.extend(["-f", "2"])

    # 执行筛选
    print(command)
    out = subprocess.run(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )  # 输出是sam

    # 可选在sam中排除chrM，将out转换为没有chrM的模式
    if rmchrM:
        out = subprocess.run(
            ["grep", "-v", "chrM"],
            input=out.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

    # 将sam转为bam
    out = subprocess.run(
        ["samtools", "view", "-h", "-b"], input=out.stdout, capture_output=True
    )

    with open(output_bam, "wb") as f:
        f.write(out.stdout)
    return para_prefix, output_bam


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


def sort_and_index_bam(bam_file, threads=5):
    basename = os.path.splitext(bam_file)[0]
    bai_file = basename + "_sorted.bai"
    sorted_bam = basename + "_sorted.bam"

    # sort bam_file
    command_sort = [
        "samtools",
        "sort",
        "-O",
        "BAM",
        "-o",
        sorted_bam,
        "-@",
        str(threads),
        bam_file,
    ]
    subprocess.run(command_sort)
    os.remove(bam_file)

    # index bam_file
    command_index = [
        "samtools",
        "index",
        sorted_bam,
        "-@",
        str(threads),
        "-o",
        bai_file,
    ]
    subprocess.run(command_index)

    return sorted_bam, bai_file


def main():
    args = parse_arguments()
    filter_para, filter_bam_file = do_bam_filter(
        bam_file=args.bam_file,
        prefix=args.prefix,
        rmdup=args.rmdup,
        mpQ_threhold=args.mpQ,
        rmUnpro_pair=args.rmunpro_pair,
        rmUnmap=args.rmunmap,
        rmchrM=args.rmchrM,
        sample_name=args.sample_name,
    )
    get_filter_qc_report(
        input_bam=filter_bam_file,
        prefix=args.prefix,
        sample_name=args.sample_name,
        filter_para=filter_para,
    )
    sorted_bam, bai_file = sort_and_index_bam(filter_bam_file, threads=args.threads)
    pass


# rmdup

# discordant

# 提取chrm

# 得到chrm比例


# 转换到bam


if __name__ == "__main__":
    main()


"""
# remove duplicates
java \
    -jar /share/home/lushaorong/apps/picard.jar MarkDuplicates \
    -I ./bowtie2/${sample_name}.bowtie2.align.sortQ30.rmChrM.bam \
    -O ./bowtie2/${sample_name}.bowtie2.align.sortQ30.rmChrMandDup.bam \
    -M ./bowtie2/${sample_name}.rmdup.matrix \
    -REMOVE_DUPLICATES true 
    # here duplicates are removed, but this may display some problems when library generation
"""
