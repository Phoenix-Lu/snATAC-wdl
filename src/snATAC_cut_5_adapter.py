# Author: Phoenix_Lu

import sys
import os
import argparse
import subprocess


def parse_arguments():
    parser = argparse.ArgumentParser(prog="Phoenix Cut 5 adaptor.")
    parser.add_argument('--r1', required=True)
    parser.add_argument('--r2', default=None)
    parser.add_argument('--mode', required=True, choices=['R1', 'R2', 'R1R2'])
    parser.add_argument('--len', required=True, type=int, nargs = '*')
    parser.add_argument('--prefix', default='_5_cutted', type=str)
    parser.add_argument('--sample_name',required=True)
    args = parser.parse_args()
    return args

def do_cut(fastq, cut_len, out_dir, prefix, sample_name):
    file_basename = os.path.basename(fastq)
    
    file_basename_without_extention = file_basename.replace('.gz', '').replace('.fq', '')
    file_out_path_without_gz = os.path.join(out_dir, sample_name + prefix)
    file_out_path = os.path.join(file_out_path_without_gz + '.fq'+ '.gz')

    command = [
        'cutadapt',
        '--cut', str(cut_len),
        '-j', '0',
        '-Z',
        '--output', str(file_out_path),
        str(fastq)
    ]

    stdout_path = file_out_path_without_gz + '.stdout'
    stderr_path = file_out_path_without_gz + '.stderr'
    out = subprocess.run(command, capture_output = True, encoding='utf-8')
    with open(stdout_path,'w') as f:
        f.write(out.stdout)
    with open(stderr_path,'w') as f:
        f.write(out.stderr)
    return file_out_path


# common usage:pinpython3 ~/work/snATAC-dev/snATAC_trim_5_adapter.py 
# --r1 ./P-single-Q5_raw_1.fq --r2 ./P-single-Q5_raw_2.fq.gz --mode 'R1R2' --len 10 10
# length can be 2 Ints
def main():
    args = parse_arguments()
    print(("#"*20) + '\n' + "do 5' cutting" + '\n' + ("#"*20))

    if (args.mode in ['R1', 'R2'] and len(args.len) == 2):
        raise ValueError("提供了大于1个cut_len")

    outdir = os.path.dirname(args.r1) + '/5_cut'
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        
    print(args.len)

    if 'R1' in args.mode:
        cut_r1_path  = do_cut(args.r1, args.len[0], outdir, args.prefix, args.sample_name)
        print(f"# r1 cutted {args.len[0]} bp at {cut_r1_path} for 5' cutting")

    if 'R2' in args.mode:
        cut_r2_path = do_cut(args.r2, args.len[-1], outdir, args.prefix, args.sample_name)
        print(f"# r2 cutted {args.len[-1]} bp at {cut_r2_path} for 5' cutting")


    print(("@"*20) + '\n' + "5' cutting Done" + '\n' + ("@"*20))


if __name__ == '__main__':
    main()
