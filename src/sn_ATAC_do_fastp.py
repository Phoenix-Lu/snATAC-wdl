# Author: Phoenix_Lu

import sys
import os
import argparse
import subprocess


def parse_arguments():
    parser = argparse.ArgumentParser(prog="Phoenix Do Fastp")
    parser.add_argument('--r1', required=True)
    parser.add_argument('--r2', default=None)
    parser.add_argument('--do_trim', action='store_true')
    parser.add_argument('--sample_name', required=True)
    parser.add_argument('--prefix', default='_trimed.fq', type=str)
    args = parser.parse_args()
    return args

def do_fastp(read1, read2, do_trim, is_pair_end, sample_name, prefix):
    print(do_trim)
    if do_trim:
        out_para = ['--out1', sample_name + '_r1' + prefix] + (['--out2', sample_name + '_r2' + prefix] if is_pair_end else [])
    else:
        out_para = []
        
    if is_pair_end:
        in_para = ['--in1', read1, '--in2', read2]
    else:
        in_para = ['--in1', read1]

    if do_trim == False :
        out_para.append('--disable_adapter_trimming')
    report_prefix = sample_name + ('_with_trim_' if (do_trim) else '_no_trim_') + 'fastp_report'
        

    command = [
        'fastp', 
        '--html', report_prefix +'.html',
        '--json', report_prefix +'.json',
        '--report_title', report_prefix,
        '--compression', '4',
        ] + in_para + out_para
    
    print(' '.join(command))

    stdout_path = report_prefix + '.stdout'
    stderr_path = report_prefix + '.stderr'
    out = subprocess.run(command, capture_output = True, encoding='utf-8')
    with open(stdout_path,'w') as f:
        f.write(out.stdout)
    with open(stderr_path,'w') as f:
        f.write(out.stderr)
    return stdout_path


# common usage:pinpython3 ~/work/snATAC-dev/snATAC_trim_5_adapter.py 
# --r1 ./P-single-Q5_raw_1.fq --r2 ./P-single-Q5_raw_2.fq.gz --mode 'R1R2' --len 10 10
# length can be 2 Ints
def main():
    args = parse_arguments()
    print(("#"*20) + '\n' + "do fastq" + '\n' + ("#"*20))
    is_pair_end = (True if args.r2 else False)

    std_out = do_fastp(read1 = args.r1,
              read2 = args.r2,
              do_trim = args.do_trim,
              is_pair_end = is_pair_end,
              sample_name = args.sample_name,
              prefix = args.prefix)


    print(("@"*20) + '\n' + "do fastq Done" + '\n' + ("@"*20))


if __name__ == '__main__':
    main()
