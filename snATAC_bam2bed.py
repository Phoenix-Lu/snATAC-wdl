# Author: Phoenix_Lu

import sys
import os
import argparse
import subprocess


def parse_arguments():
    parser = argparse.ArgumentParser(prog="Phoenix Bowtie2 Align adaptor.")
    parser.add_argument("--pair_ends", action="store_true")
    parser.add_argument("--do_atac_shift", default=False, action="store_true")
    parser.add_argument("--bam_file", required=True)
    parser.add_argument("--black_list", type=str, default=None)
    parser.add_argument("--threads", default=5)
    parser.add_argument("--prefix", default="_bowtie2Mapped", type=str)
    parser.add_argument("--sample_name", required=True)
    args = parser.parse_args()
    return args


def bam2bed(bam_file, black_list=None, do_atac_shift=False):
    tmp_bed = "tmp_.bed"
    basename = os.path.basename(bam_file)
    basename = os.path.splitext(basename)[0]
    bed_file = basename + ".bed"
    command = ["bedtools", "bamtobed", "-i", bam_file, ">", tmp_bed]

    print(command)
    subprocess.run(" ".join(command), shell=True, capture_output=True)

    if black_list:
        command = [
            "bedtools",
            "intersect",
            "-v",
            "-a",
            tmp_bed,
            black_list,
            ">",
            bed_file,
        ]
        subprocess.run(" ".join(command), shell=True, capture_output=True)

#TODO debug shift
    if do_atac_shift:
        command = ["cat", tmp_bed, "|", r"""awk -F '\t' 'BEGIN {OFS = FS}{ \
  if ($6 == "+") {$2 = $2 + 4} \
  else if ($6 == "-") {$3 = $3 - 5} \
  print $0}'""",">", bed_file]
        print(' '.join(command))
        # command = (
        #     f"""cat {tmp_bed}"""
        #     + r""" | awk -F '\t'"""
        #     + """'BEGIN {OFS = FS}{if ($6 == \'+\') {$2 = $2 + 4} else if ($6 == \'-\') {$3 = $3 - 5} print $0}' > """
        #     + f"""{bed_file}"""
        # )
        print(command)
        subprocess.run(command, shell=True)

    else:
        subprocess.run(f"mv {tmp_bed} {bed_file}", shell=True)

    return bed_file


def main():
    args = parse_arguments()
    bam2bed(bam_file=args.bam_file, do_atac_shift=True)


if __name__ == "__main__":
    main()
