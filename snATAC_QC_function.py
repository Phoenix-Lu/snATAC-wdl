import subprocess
import os


def get_bam_qc_report(input_bam, prefix, sample_name, filter_para=None):
    qc_report = sample_name + prefix + ".filtersummary"

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
