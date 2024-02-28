# Author: Phoenix_Lu

import sys
import os
import argparse
import subprocess

# # mapping
# bowtie2 \
#     -x ${ref_path_bowtie2} \
#     --very-sensitive \
#     -1 ./cutadapt.cut32/${sample_name}_f1.trimed.fq \
#     -2 ./cutadapt.cut32/${sample_name}_r2.trimed.fq \
#     -p 20 \
#     -X 2000 \
#     --no-mixed \
#     --no-discordant \
#     -S ./bowtie2/${sample_name}.bowtie2.align.sam \
#     2> ./${sample_name}.bowtie2.align.summary


def parse_arguments():
    parser = argparse.ArgumentParser(prog="Phoenix Bowtie2 Align adaptor.")
    parser.add_argument("--r1", required=True)
    parser.add_argument("--r2", type=str, default=None)
    parser.add_argument("--reference_path", required=True)
    parser.add_argument("--mode", default="very-sensitive", choices=["very-sensitive"])
    parser.add_argument("--threads", default=5)
    parser.add_argument(
        "--prefix", default="_bowtie2Mapped", type=str
    )  # 对于type=str, 如果输入参数则会发生改变
    parser.add_argument(
        "--no_discordant", default=False, action="store_true"
    )  # 对于此类TF，只能默认default=False，无法指定--para False
    parser.add_argument("--sample_name", required=True)
    parser.add_argument("--pair_ends", action="store_true")
    parser.add_argument("--preserve_sam", default=False, action="store_true")
    args = parser.parse_args()
    return args


# PASSED
def do_bowtie2_align_PE(
    r1,
    r2,
    mode,
    reference_path,
    sample_name,
    threads,
    no_discordant,
    prefix="_bowtie2Mapped",
):

    sam_file_path = sample_name + prefix + ".sam"
    command = [
        "bowtie2",
        "-x",
        str(reference_path),
        "-1",
        str(r1),
        "-2",
        str(r2),
        "-p",
        str(threads),
        "--" + str(mode),
        "-X",
        str(2000),
        "-S",
        str(sam_file_path),
    ]
    if no_discordant:
        command.extend(["--no_discordant"])
    print(" ".join(command))
    stderr_path = sample_name + prefix + ".alignsummary"
    out = subprocess.run(command, capture_output=True, encoding="utf-8")
    with open(stderr_path, "a") as f:
        print(f"{'*' * 10} {sample_name} Bowtie2 report {'*' * 10}\n", file=f)
        print(
            f"Parameters:\n \
              ref_path:{reference_path}\n \
                r1:{r1}\n \
                    r2:{r2}\n \
                        mode:--{mode}\n \
                            maxlen:2000\n \
                                samoutput:{sam_file_path}\n \
                                    {'--no_discordant' if no_discordant else '' }",
            file=f,
        )
        f.write(out.stderr)
    return sam_file_path


# 用samtools stat产生报告，并且获取我需要的值


def trans_sam_to_sorted_bam(sam_file, threads, preserve_sam):
    basename = os.path.splitext(sam_file)[0]
    tmp_bam = str(basename) + "_tmp" + ".bam"
    sorted_bam = str(basename) + ".bam"
    command_trans = [
        "samtools",
        "view",
        sam_file,
        "-O",
        "BAM",
        "-o",
        tmp_bam,
        "-h",
        "-@",
        str(threads),
    ]
    subprocess.run(command_trans)

    command_sort = [
        "samtools",
        "sort",
        "-O",
        "BAM",
        "-o",
        sorted_bam,
        "-@",
        str(threads),
        tmp_bam,
    ]
    subprocess.run(command_sort)
    os.remove(tmp_bam)
    if not preserve_sam:
        os.remove(sam_file)
    return sorted_bam


def picard_markdup(
    input_bam, remove_dup=False, picard_path=r"/share/home/lushaorong/apps/picard.jar"
):
    basename = os.path.splitext(input_bam)[0]
    rmdup_bam = basename + ("_markdup" if remove_dup == False else "_rmdup") + ".bam"
    command = [
        "java",
        "-jar",
        picard_path,
        "MarkDuplicates",
        "-I",
        input_bam,
        "-O",
        rmdup_bam,
        "-M",
        basename + "_rmdpup_matrix",
    ]
    if remove_dup:
        command.extend(["-REMOVE_DUPLICATES", "true"])
    subprocess.run(command)
    return rmdup_bam


def get_align_qc_report(input_bam, prefix, sample_name, filter_para=None):
    align_qc_report = sample_name + prefix + ".full_alignsummary"
    # 打印stat 报告
    with open(align_qc_report, "a") as f:
        print(
            f"{'*' * 20} {input_bam} full mapping QC Report by samtools stat {'*' * 20}",
            file=f,
        )
        print(filter_para, file=f)
    command = f"samtools stat {input_bam} |grep ^SN | cut -f 2- >> {align_qc_report}"
    subprocess.run(command, shell=True)

    # 打印flagstat报告
    with open(align_qc_report, "a") as f:
        print("\n\n", file=f)
        print(
            f"{'*' * 20} {input_bam} full mapping QC Report by samtools flagstat {'*' * 20}",
            file=f,
        )
    command = f"samtools flagstat {input_bam} >> {align_qc_report}"
    subprocess.run(command, shell=True)
    return align_qc_report


def get_align_qc_report_summary(full_qc_report, prefix, sample_name):
    align_qc_report = sample_name + prefix + ".alignsummary"


def main():

    args = parse_arguments()
    # 如果没有输入r2但是指定双端
    is_pair_ends = args.pair_ends
    if args.pair_ends and (not args.r2):
        raise TypeError("--pair_ends defined but r2 not provided")
    if args.r1 and args.r2:
        is_pair_ends = True
    # align
    if is_pair_ends == True:
        print(("#" * 20) + "\n" + "do PE Align" + "\n" + ("#" * 20) + "\n")
        out_sam_file = do_bowtie2_align_PE(
            r1=args.r1,
            r2=args.r2,
            reference_path=args.reference_path,
            mode=args.mode,
            threads=args.threads,
            no_discordant=args.no_discordant,
            sample_name=args.sample_name,
            prefix=args.prefix,
        )

    else:
        # Do_bowtie2_align_SE
        pass

    orgin_bam = trans_sam_to_sorted_bam(
        sam_file=out_sam_file, threads=args.threads, preserve_sam=args.preserve_sam
    )
    rmdup_bam = picard_markdup(input_bam=orgin_bam)
    full_align_qc_report = get_align_qc_report(
        input_bam=rmdup_bam, prefix=args.prefix, sample_name=args.sample_name
    )
    # rmdup

    # discordant

    # 提取chrm

    # 得到chrm比例

    # 转换到bam
    print(
        ("@" * 20)
        + "\n"
        + f"Aligning Done, alignfile at {out_sam_file}"
        + "\n"
        + ("@" * 20)
    )


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
