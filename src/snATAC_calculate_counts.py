# Author: Phoenix_Lu

import sys
import os
import argparse
import subprocess
import pandas
from io import StringIO


def parse_arguments():
    parser = argparse.ArgumentParser(prog="Phoenix caculate cell counts")
    parser.add_argument('--queue_bam', required=True, nargs= '+', help="BED is also ok but i recoomend")
    parser.add_argument('--peak_reference_bed', required=True, help="can be pooled called_peaks or general reference")
    parser.add_argument('--experiment_name', required=True)
    args = parser.parse_args()
    return args

# 似乎使用bam文件更为靠谱， 如果是bed，若单个read在范围内会被计算，所有insert或者覆盖度相关的都从bam开始算。
# 结论：用bam而不是bed。bed只作为区域指定的文件。
# coverageBed的结果怎么解读？
    # 1) depth ：最需要的
    # 2) bases at depth ： 
    # 3) size of A ： reference
    # 4) % of A at depth

# callpeak可以不严格

# TEST_PASSED
def count_one_cell(queue_bam, reference_bed, cell_id):
    # A是reference， B是queue
    command = ['coverageBed', '-a', reference_bed, '-b', queue_bam, "|", "cut", "-f", "4,5", ">", f"{cell_id}.count"] # -b in -a
    print(' '.join(command))
    subprocess.run(' '.join(command), shell=True)
    return str(f"{cell_id}.count")
    # #让intersect看起来像个文件
    # intersect = StringIO(intersect)
    # pd = pandas.read_csv(intersect, sep='\t',header=None)
    # cell_peak_index = pd[3]
    # cell_peak_depth = pd[4]
    # return {"cell_id": cell_id, "peak_index": cell_peak_index, "cell_peak_depth" : cell_peak_depth}

# TEST_PASSED
def gen_count_matrix(sc_count:list, ref_bed, experiment_name):
    # output csc_matrix format
    ## rows are peaks, columns are cells
    import subprocess
    header = '%%MatrixMarket matrix coordinate integer general\n'
    spacer = '%\n'

    mtx = open(f'{experiment_name}_count.mtx', 'w')
    mtx.write(header)
    mtx.write(spacer)

    rownames = open(f'{experiment_name}_count.rows', 'w')
    colnames = open(f'{experiment_name}_count.cols', 'w')

    # 把peak的位置信息写入.rows中
    with open(ref_bed) as peak:
        for line in peak:
            pos = line.split('\t')[:3]
            rownames.write('\t'.join(pos))
            rownames.write('\n')

    total_non_zeroes = 0
    # id, iter_thing in enumerate()
    # c is cell
    for c, j in enumerate(sc_count):
        # pl is sample
        pl = experiment_name
        # pl = "sample"
        # fn is basename
        fn = j
        colnames.write(pl + '_' + fn[:-6] + '\n')
        # p is peaks
        with open(j) as fh:
            for p, v in enumerate(fh):
                pid, count = v.strip().split('\t')
                coln = c + 1
                if count != '0':
                    total_non_zeroes += 1
                    rown = p + 1
                    mtx.write('%s %s %s\n' %(rown, coln, count))

    total_c = c + 1
    total_p = p + 1
    print(f"total cells: {total_c}\ntotal peaks: {total_p}\ntotal_non_zeroes: {total_non_zeroes}")


    # add the information of numbers of rows, columns and non-zeroes
    # to the third (the "3i" command of sed) line of the mtx file
    # subprocess.run("touch count_matrix_over_aggregate.mtx", shell=True)
    mtx.close()
    rownames.close()
    colnames.close()
    subprocess.run(f"sed -i '3i%s %s %s' {experiment_name}_count.mtx" %
                (total_p, total_c, total_non_zeroes),
                shell=True)



# queue_bam = "MEF_sub_filtered_sorted.bam"
# reference_bed = "MEF_sub_filtered_sorted_summits.bed"
# count_one_cell(queue_bam, reference_bed, "cell_1")
# count_one_cell(queue_bam, reference_bed, "cell_2")
# sc_count = ["cell_1.count", "cell_2.count"]
# gen_count_matrix(sc_count, reference_bed)


def main():
    args = parse_arguments()
    cell_count = []
    for cell in args.queue_bam:
        basename = (os.path.basename(cell)).split('.')[0]
        cell_count.append(count_one_cell(queue_bam=cell, reference_bed=args.peak_reference_bed,cell_id = basename))
    print("sc Count generated")
    gen_count_matrix(sc_count=cell_count, ref_bed=args.peak_reference_bed,experiment_name= args.experiment_name)



if __name__ == '__main__':
    main()


