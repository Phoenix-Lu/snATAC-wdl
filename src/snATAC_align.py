# Author: Phoenix_Lu

import sys
import os
import argparse
import subprocess

# TESTPASSED
# 1.标准带r1r2已完成测试
# 2.preserve_sam ok
# 3.preserve_raw_bam ok


# FIXME align使用的reference需要与后续bam2bed一致
def parse_arguments():
    parser = argparse.ArgumentParser(prog="Phoenix Bowtie2 Align adaptor.")
    parser.add_argument("--r1", required=True)
    parser.add_argument("--r2", type=str, default=None)
    parser.add_argument("--sample_name", required=True)
    parser.add_argument("--picard_path",default = r"/share/home/lushaorong/apps/picard.jar")
    parser.add_argument("--reference_path", required=True)


    parser.add_argument("--no_discordant", default=False, action="store_true") 
    # 对于此类TF，只能默认default=False，无法指定--para False
    parser.add_argument("--mode", default="very-sensitive", choices=["very-sensitive", "sensitive", "fast", "very-fast"])
    parser.add_argument("--threads", default=20)
    parser.add_argument("--prefix", default="_bowtie2Mapped_raw", type=str) 

    # 对于type=str, 如果输入参数则会发生改变
    parser.add_argument("--pair_ends", action="store_true")
    parser.add_argument("--preserve_sam", default=False, action="store_true")
    parser.add_argument("--preserve_raw_bam", default=False, action="store_true")

    parser.add_argument("--remove_dup", default=False, action="store_true")
    args = parser.parse_args()
    return args


# PASSED
def do_bowtie2_align_PE(r1, r2, mode, reference_path, sample_name, threads, no_discordant, prefix):
    sam_file_path = sample_name + prefix + ".sam"
    command = [
        "bowtie2",
        "-x", str(reference_path),
        "-1", str(r1),
        "-2", str(r2),
        "-p", str(threads),
        "--" + str(mode),
        "-X", str(2000),
        "-S", str(sam_file_path)]
    if no_discordant:
        command.extend(["--no_discordant"])

    print(" ".join(command))
    out = subprocess.run(command, capture_output=True, encoding="utf-8")

    # stderr保存了align
    stderr_path = sample_name + prefix + ".alignsummary"
    with open(stderr_path, "a") as f:
        print(f"{'*' * 10} {sample_name} Bowtie2 report {'*' * 10}\n", file=f)
        # 打印所有参数
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


def picard_markdup(input_bam, remove_dup=False, picard_path=r"/share/home/lushaorong/apps/picard.jar"):
    basename = os.path.splitext(input_bam)[0]
    rmdup_bam = basename + ("_markdup" if remove_dup == False else "_rmdup") + ".bam"
    command = [
        "java", "-jar", picard_path,
        "MarkDuplicates",
        "-I", input_bam,
        "-O", rmdup_bam,
        "-M", basename + ".picardDupsummary"]
    if remove_dup:
        command.extend(["-REMOVE_DUPLICATES", "true"])
    subprocess.run(command)
    return rmdup_bam, basename + ".picardDupsummary"



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
        # TODO Do_bowtie2_align_SE
        pass
    
    # sam转bam
    orgin_bam = trans_sam_to_sorted_bam(sam_file=out_sam_file,
                                         threads=args.threads,
                                           preserve_sam=args.preserve_sam)
    
    # picard标重
    rmdup_bam, picardDupsummary= picard_markdup(input_bam=orgin_bam,remove_dup = args.remove_dup)

    # TODO 得到alignreport，质检需要
    # full_align_qc_report = get_align_qc_report(input_bam=rmdup_bam,
    #                                             prefix=args.prefix,
    #                                               sample_name=args.sample_name)
    if not args.preserve_raw_bam:
        os.remove(orgin_bam)

    print("\n" + ("@" * 20) + "\n"+
           f"Aligning Done,\n"+
            (f"alignRawBam at {orgin_bam}\n" if (args.preserve_raw_bam) else "alignRawBam is delected\n")+
           f"marked_RawBam at {rmdup_bam}\n" +
           f"picard duplicate metric at {picardDupsummary}\n"+
             ("@" * 20))


if __name__ == "__main__":
    main()


