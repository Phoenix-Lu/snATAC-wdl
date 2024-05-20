# Author: Phoenix_Lu

import sys
import os
import argparse
import subprocess


# TEST PASSED
# 1.r1 r2均提供，do or not trim 可通过
# 2.只提供r1，do or not trim均通过，但是不会trim
# 返回的是.fq.gz

# 注意:fastqc的trim主要基于overlap判断adapter
#FIXME SE的情况下，adapter检查不佳，发现似乎不会去主动删除adapter
def parse_arguments():
    parser = argparse.ArgumentParser(prog="Phoenix Do Fastp")
    parser.add_argument('--r1', required=True)
    parser.add_argument('--r2', default=None)
    parser.add_argument('--do_trim', action='store_true')
    parser.add_argument('--sample_name', required=True)
    parser.add_argument('--prefix', default='_trimed.fq.gz', type=str)
    args = parser.parse_args()
    return args

def do_fastp(read1, read2, do_trim, is_pair_end, sample_name, prefix):  
    print(f"do_trim? {do_trim}")

    if do_trim:
        out_para = ['--out1', sample_name + '_r1' + prefix] + (['--out2', sample_name + '_r2' + prefix] if is_pair_end else [])
    else:
        out_para = []
        out_para.append('--disable_adapter_trimming')
        
    if is_pair_end:
        in_para = ['--in1', read1, '--in2', read2]
    else:
        in_para = ['--in1', read1]
        
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
    stderr_path = report_prefix + '.fastpsummary'
    out = subprocess.run(command, capture_output = True, encoding='utf-8')
    # 若有信息再写入
    if out.stdout:
        with open(stdout_path,'w') as f:
            f.write(out.stdout)
    if out.stderr:
        with open(stderr_path,'w') as f:
            f.write(out.stderr)
    return stderr_path


def main():
    args = parse_arguments()
    print(("#"*20) + '\n' + "do fastq" + '\n' + ("#"*20) + '\n')
    is_pair_end = (True if args.r2 else False)

    stderr_path = do_fastp(read1 = args.r1,
              read2 = args.r2, 
              do_trim = args.do_trim,
              is_pair_end = is_pair_end,
              sample_name = args.sample_name,
              prefix = args.prefix)


    print('\n' + ("@"*20) + '\n' + "do fastq Done" + '\n' + ("@"*20) + '\n')


if __name__ == '__main__':
    main()
