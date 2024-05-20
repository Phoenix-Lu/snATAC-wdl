# Author: Phoenix_Lu

import sys
import os
import argparse
import subprocess
import pandas

# TESTPASSED
# do_atac_shift和black_list 交叉验证通过
def parse_arguments():
    parser = argparse.ArgumentParser(prog="Phoenix Bam2Bed.")
    # required
    parser.add_argument("--bam_file", required= True)

    # bw_needed
    parser.add_argument("--chrom_sizes")
    parser.add_argument("--genome_bed")

    # process_arg
    parser.add_argument("--do_atac_shift", default=False, action="store_true")
    parser.add_argument("--black_list", type=str, default=None)
    # eff_arg
    parser.add_argument("--threads", default=5)


    args = parser.parse_args()
    return args


def bam2bed(bam_file, black_list=None, do_atac_shift=False,
               chrom_sizes = None, genome_bed = None):
    tmp_bed = "tmp_.bed"
    basename = os.path.basename(bam_file)
    basename = os.path.splitext(basename)[0]
    shift_tmp_bed = 'shift_bed.tmp'
    bed_file = basename + ".bed"
    command = ["bedtools", "bamtobed", "-i", bam_file, ">", tmp_bed]
    print("***bam ---> bed***")
    print(" ".join(command))
    subprocess.run(" ".join(command), shell=True, capture_output=True)
    # 以上得到tmp.bed
    
    # shift功能通过测试
    # tmp.bed -> shift_tmp_bed
    if do_atac_shift:
        print("***DO ATAC SHIFT***")
        tmp_bed_pd = pandas.read_csv(tmp_bed, sep="\t", header=None)
        # 正链+4，负链-5
        tmp_bed_pd[1] = tmp_bed_pd.apply(lambda row: row[1] + 4 if row[5] == '+' else row[1], axis=1)
        tmp_bed_pd[2] = tmp_bed_pd.apply(lambda row: row[2] - 5 if row[5] == '-' else row[2], axis=1)
        tmp_bed_pd.to_csv(shift_tmp_bed, sep='\t', header=False, index=False)
        os.remove(tmp_bed)

    else:
        print("***NO ATAC SHIFT***")
        print(f"mv {tmp_bed}  {shift_tmp_bed}")
        subprocess.run(f"mv {tmp_bed}  {shift_tmp_bed}", shell=True)

    if black_list:
        command = [
            "bedtools",
            "intersect",
            "-v",
            "-a",
            shift_tmp_bed,
            "-b",
            black_list,
            ">",
            bed_file,
        ]
        print("***DO blacklist filter***")
        print(' '.join(command))
        subprocess.run(" ".join(command), shell=True, capture_output=True)
        os.remove(shift_tmp_bed)
    else:
        print("***NO blacklist filter***")
        command  = f"mv {shift_tmp_bed} {bed_file}"
        print(command)
        subprocess.run(command, shell=True)


    if chrom_sizes and genome_bed:
        bw_file = bed2bigwig(bed_files=bed_file, genome_bed=genome_bed, chrom_sizes=chrom_sizes)
    else:
        print("lack file to generate .bw file, skipped")
        print(f"""lack {"chrom_sizes" if not chrom_sizes else ""} {"genome_bed" if not genome_bed else ""}""")
        bw_file = None
    return bed_file, bw_file



#TODO bed2bigwig 目测ok，需要串起来检查,目前缺文件，记得找时间重新搭建reference

def bed2bigwig(bed_files, genome_bed, chrom_sizes):
    # genome_file is .fai
    basename = os.path.basename(bed_files)
    basename = os.path.splitext(basename)[0]
    bed_graph_file = basename + '.bedgraph'
    bigwig_file = basename+'.bw'
    #  1052  2024-03-04 15:01:41 sort -k1,1V -k2,2n -k3,3n filter_all_filtered_sorted.bed > filter_all_filtered_sorted.sorted.bed
    
    commmand = ['sort -k 1,1', bed_files,'>', 'tmp_bed', ";",
        'bedtools', 'genomecov', '-i', 'tmp_bed', '-g', chrom_sizes, '-bg','>', bed_graph_file, ';',
                 'bedGraphToBigWig', bed_graph_file, chrom_sizes, bigwig_file]
    print(' '.join(commmand))
    subprocess.run(' '.join(commmand), shell=True)
    os.remove('tmp_bed')
    
def main():
    args = parse_arguments()
    bam2bed(bam_file=args.bam_file, do_atac_shift=args.do_atac_shift, black_list=args.black_list,
            chrom_sizes= args.chrom_sizes, genome_bed= args.genome_bed)


if __name__ == "__main__":
    main()
