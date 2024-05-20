# Author: Phoenix_Lu

import sys
import os
import argparse
import subprocess


def parse_arguments():
    parser = argparse.ArgumentParser(prog="Phoenix Do Fastp")
    parser.add_argument('--bam_files', required=True, nargs='+') 
    # nargs = '+' 可以帮助传入多个参数，传入参数时，通过 --arg A B C --arg_next
    parser.add_argument('--experiment_name', required=True)
    parser.add_argument('--prefix', default='_aggregate.bam', type=str)
    args = parser.parse_args()
    return args

def aggregate_bam(bam_files:list, prefix, experiment):
    aggregate_bam = experiment + prefix
    command = ['samtools', 'merge',
               '-o', aggregate_bam]+bam_files
    print(command)
    subprocess.run(command)
    return aggregate_bam

def main():
    args = parse_arguments()
    print(("#"*20) + '\n' + f"do merge for {args.experiment_name}" + '\n' + ("#"*20))
    aggregated_bam_file = aggregate_bam(bam_files=args.bam_files, prefix=args.prefix, experiment = args.experiment_name)


    print(("@"*20) + '\n' + "do aggregate Done" + '\n' + ("@"*20))
    return aggregate_bam

if __name__ == '__main__':
    main()
