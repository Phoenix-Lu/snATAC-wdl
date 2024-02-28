# Author: Phoenix_Lu

import sys
import os
import argparse
import subprocess

def parse_arguments():
    parser = argparse.ArgumentParser(prog="Phoenix Bam filter.")
    parser.add_argument('--bam_file', required=True)
    parser.add_argument('--mpQ', default = '30')
    parser.add_argument('--rmunpair', rdafaultequired = True)
    parser.add_argument('--rmunmap', dafault=True)
    parser.add_argument('--rmdup', dafault=True)
    parser.add_argument('--threads', default = 5)
    parser.add_argument('--prefix', default='_filtered', type=str)
    parser.add_argument('--sample_name', required=True)
    args = parser.parse_args()
    return args

def do_bam_filter(bam_file, prefix):
    basename = os.path.splitext(bam_file)[0]
    filtered_bam = basename + prefix + '.bam'
    





def do_bowtie2_align_PE(r1 , r2, prefix, mode, reference_path,
                           sample_name, threads, no_discordant):
    sam_file_path = sample_name + prefix + '.sam'
    command = [
        'bowtie2',
        '-x', str(reference_path),
        '-1', str(r1),
        '-2', str(r2),
        '-p', str(threads), 
        '--' + str(mode),
        '-X', str(2000),
        '-S', str(sam_file_path)]
    if (no_discordant):
        command.extend(['--no_discordant'])        
    print(' '.join(command))
    stderr_path = sample_name + prefix + '.alignsummary'

    out = subprocess.run(command, capture_output = True, encoding='utf-8')
    with open(stderr_path,'a') as f:
        print(f'{'*' * 10} {sample_name} Bowtie2 report {'*' * 10}\n', file = f)
        print(f'Parameters:\n \
              ref_path:{reference_path}\n \
                r1:{r1}\n \
                    r2:{r2}\n \
                        mode:--{mode}\n \
                            maxlen:2000\n \
                                samoutput:{sam_file_path}\n \
                                    {'--no_discordant' if no_discordant else ''}', file=f)
        f.write(out.stderr)
    return sam_file_path


# 用samtools stat产生报告，并且获取我需要的值

def trans_sam_to_sorted_bam(sam_file, threads, preserve_sam):
    basename = os.path.splitext(sam_file)[0]
    tmp_bam = str(basename) + '_tmp' + '.bam'
    sorted_bam = str(basename) + '.bam'
    command_trans = [
        'samtools', 'view',
        sam_file,
        '-O', 'BAM',
        '-o', tmp_bam,
        '-h',
        '-@', str(threads)]
    subprocess.run(command_trans)

    command_sort = [
        'samtools', 'sort',
        '-O', 'BAM',
        '-o', sorted_bam,
        '-@', str(threads),
        tmp_bam]
    subprocess.run(command_sort)
    os.remove(tmp_bam)
    if not preserve_sam:
        os.remove(sam_file)
    return sorted_bam

def picard_markdup(input_bam, remove_dup = False, picard_path = r'/share/home/lushaorong/apps/picard.jar'):
    basename = os.path.splitext(input_bam)[0]
    rmdup_bam = basename + ('_markdup' if remove_dup ==False else '_rmdup') + '.bam'
    command = ['java', '-jar',
               picard_path,
               'MarkDuplicates',
               '-I', input_bam,
               '-O', rmdup_bam,
               '-M', basename + '_rmdpup_matrix'
               ]
    if remove_dup:
        command.extend(['-REMOVE_DUPLICATES', 'true'])
    subprocess.run(command)
    return rmdup_bam
    
def get_align_qc_report(input_bam, prefix, sample_name):
        align_qc_report = sample_name + prefix + '.alignsummary'
        
        with open(align_qc_report, 'a') as f:
            print(f'{'*' * 10} {input_bam} full mapping QC Report {'*' * 10}', file = f)
        command = f'samtools stat {input_bam} |grep ^SN | cut -f 2- >> {align_qc_report}'
        subprocess.run(command, shell=True)
        return align_qc_report




def main():
    
    args = parse_arguments()
    # 如果没有输入r2但是指定双端
    is_pair_ends = args.pair_ends
    if args.pair_ends and (not args.r2):
        raise TypeError('--pair_ends defined but r2 not provided')
    if args.r1 and args.r2:
        is_pair_ends = True
# align
    if is_pair_ends == True:
        print(("#"*20) + '\n' + "do PE Align" + '\n' + ("#"*20) + '\n')
        out_sam_file = do_bowtie2_align_PE(r1 = args.r1, r2 = args.r2, reference_path = args.reference_path,
                                            mode = args.mode, threads = args.threads,
                                            no_discordant= args.no_discordant,
                                              sample_name = args.sample_name, 
                                              prefix = args.prefix)
                                              
    else:
        # Do_bowtie2_align_SE
        pass

    orgin_bam = trans_sam_to_sorted_bam(sam_file = out_sam_file, threads = args.threads, preserve_sam=args.preserve_sam)
    rmdup_bam = picard_markdup(input_bam=orgin_bam)
# rmdup

# discordant

# 提取chrm
    
# 得到chrm比例
    

# 转换到bam
    print(("@"*20) + '\n' + f"Aligning Done, alignfile at {out_sam_file}" + '\n' + ("@"*20))
    
    

if __name__ == '__main__':
    main()


"""
# remove duplicates
java \
    -jar /share/home/lushaorong/apps/picard.jar MarkDuplicates \
    -I ./bowtie2/${sample_name}.bowtie2.align.sortQ30.rmChrM.bam \
    -O ./bowtie2/${sample_name}.bowtie2.align.sortQ30.rmChrMandDup.bam \
    -M ./bowtie2/${sample_name}.rmdup.matrix \
    -REMOVE_DUPLICATES true 
    # here duplicates are removed, but this may display some problems when library generation
"""