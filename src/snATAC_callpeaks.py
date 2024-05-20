# Author: Phoenix_Lu

import sys
import os
import argparse
import subprocess


def parse_arguments():
    parser = argparse.ArgumentParser(prog="Phoenix macs2 call peak.")
    parser.add_argument('--bed_file', required=True)
    parser.add_argument('--genome_size', required=True, choices=['mm', 'hs'])

    parser.add_argument('--ATAC_default_mode', default=False, action = 'store_true')
    
    
    parser.add_argument('--shift', default = -75, help = "extend == -2*shift. -75 is the default arg for ATAC")
    parser.add_argument('--p_val', default = 0.01)
    
    parser.add_argument('--control', default=None)
    args = parser.parse_args()
    return args


def call_peak_macs2(bed_file, genome_size, shift, p_val, bed_control = None):
    basename = os.path.basename(bed_file).split('.')[0]
    command = [
        'macs2',
        'callpeak',
        '-t', str(bed_file),
        '-f', 'BED',
        '-g', str(genome_size),
        '-n', basename,
        '--outdir', './',
        '--shift', str(shift),
        '--nomodel',
        '--extsize', str(-2*shift),
        '--pvalue', str(p_val),
         # 保存峰轮廓和片段堆积信息与bedgraph，需要保存
        '-B', '--SPMR']
    if bed_control:
        command.extend(['-c', bed_control])

    print(' '.join(command))

    stdout_path = basename + '.stdout'
    stderr_path = basename + '.stderr'
    out = subprocess.run(command, capture_output=True, encoding='utf-8')
    if out.stdout:
        with open(stdout_path,'w') as f:
            f.write(out.stdout)
    if out.stderr:
        with open(stderr_path,'w') as f:
            f.write(out.stderr)
    # return sample_name + prefix + '_peaks.narrowPeak'


def main():
    args = parse_arguments()
    if args.ATAC_default_mode:
        args.shift = -75
        args.p_val = 0.01
        print("Using ATAC default args:")

    if args.control :
        print("Control provided") 
    else:
        print("Control not provided") 
    call_peak_macs2(bed_file = args.bed_file,
                    genome_size = args.genome_size,
                    shift = args.shift,
                    p_val = args.p_val,
                    bed_control= args.control)


if __name__ == '__main__':
    main()
