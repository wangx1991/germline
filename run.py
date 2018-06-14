import logging
import os
import re
import sys
import time
import gzip
import itertools
from difflib import SequenceMatcher
import argparse

#sys.path.append(os.path.split(os.path.realpath(__file__))[0])
#import the logging 
from pipelines.log.log import store_pipeline_logs , store_trim_logs, store_filter_logs, store_align_logs , store_germline_VC_logs
#import the trim functions
from pipelines.trim.trim_reads import trim_read_pairs
#import the align functions
from pipelines.align.align_reads import align_reads_bwa
#import the align functions
from pipelines.post_align.post_alignment import filter_alignment_samtools , identify_gs_primers


def script_information():
    print "\nApplication: pipelines of QIAseq Targeted DNA Panel\n"
    print "====================================================================="
    print "Required environment£ºpython \ bwa \ samtools \ GATK"

parser = argparse.ArgumentParser(usage = "\n\npython %(prog)s --source --sample_name --tailname --common_seq1 --common_seq2 --output --min_read_len --bwa_dir --ref_index_name --ref_fa_file --num_threads --samtools_dir --min_mapq --max_soft_clip --max_dist")
parser.add_argument("--source", help = "Path to input reads in FASTA format", type = str)
parser.add_argument("--sample_name", help = "the sample name of raw reads", type = str)
parser.add_argument("--tailname", help = "the tailname of sample raw reads", type =s tr)
parser.add_argument("--common_seq1", help = "the common seq1 of QIAseq Targeted DNA Panel", type = str, default = 'CAAAACGCAATACTGTACATT')
parser.add_argument("--common_seq2", help = "the common seq2 of QIAseq Targeted DNA Panel", type = str, default = 'ATTGGAGTCCT')
parser.add_argument("--output", help = "Path of output file", type = str)
parser.add_argument("--min_read_len", help = "the cutoff of the min read length", type = int, default = 40)
parser.add_argument("--bwa_dir", help = "the install path of bwa", type = str)
parser.add_argument("--ref_index_name", help = "the path of ref index£¬if there isn't a ref index, it will make a index in the path of ref fasta by bwa", type = str)
parser.add_argument("--ref_fa_file", help = "the path of ref fasta", type = str)
parser.add_argument("--num_threads", help = "the number of threads to align", type = int, default = 4)
parser.add_argument("--samtools_dir", help = "the install path of samtools", type = str)
parser.add_argument("--min_mapq", help = "the parameter of filter alignment_sam", type = int, default = 17)
parser.add_argument("--max_soft_clip", help = "the parameter of filter alignment_sam", type = int, default = 10)
parser.add_argument("--max_dist", help = "the parameter of filter alignment_sam", type = int, default = 2)
parser.add_argument("-v", '-version', action = 'version', version =' %(prog)s 1.0')

if "-h" in sys.argv[1:] or "--help" in sys.argv[1:]:
    script_information()
    args = parser.parse_args()
if "-v" in sys.argv[1:]:
    args = parser.parse_args()
try:
    args = parser.parse_args()
except SystemExit:
    script_information()
    parser.print_help()
    exit()

if len(sys.argv) == 1:
    script_information()
    parser.print_help()
    exit()

def main():
    #---check the outputdir
    out_dir = args.out_dir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    #----pipeline log file
    log_dir = out_dir + '/' + 'log/'
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    logger_pipeline_process, logger_pipeline_errors = store_pipeline_logs(log_dir)

    #---input
    source = args.source
    sample = args.sample_name
    tailname = args.tailname
    min_read_len = args.min_read_len
    common_seq1 = args.common_seq1
    common_seq2 = args.common_seq2
    num_threads = args.num_threads
    bwa_dir = args.bwa_dir
    ref_fa_file = args.ref_fa_file
    ref_index_name = args.ref_index_name
    samtools_dir = args.samtools_dir
    min_mapq = args.min_mapq
    max_soft_clip = args.max_soft_clip
    primers_file = args.primers_file
    ##########################################################################################
    #---trim
    ##########################################################################################
    #undetermined_dir
    undetermined_dir = out_dir + '/'+ 'undetermined'
    os.makedirs(undetermined_dir)

    log_trim_dir = undetermined_dir + '/' + 'log/'
    if not os.path.exists(log_trim_dir):
        os.makedirs(log_trim_dir)

    tailname = '_' + tailname
    read1 = source + '/' + sample + tailname + '_R1_001.fastq.gz'
    read2 = source + '/' + sample + tailname + '_R2_001.fastq.gz'
    trimmed1 = undetermined_dir + '/' + sample + '_R1_undetermined.fastq'
    trimmed2 = undetermined_dir + '/' + sample + '_R2_undetermined.fastq'
    stats_file = undetermined_dir + '/' + sample + '_basic_stats.txt'

    logger_trim_process, logger_trim_errors = store_trim_logs(log_trim_dir)
    trim_read_pairs(read1, read2, trimmed1, trimmed2, min_read_len,
                    common_seq1, common_seq2, stats_file, logger_trim_process,
                    logger_trim_errors)
                           common_seq1, common_seq2, stats_file, logger_trim_process,
                           logger_trim_errors)

    ##########################################################################################
    #---align
    ##########################################################################################
    #aligned_dir
    aligned_dir = out_dir + '/'+ 'aligned'
    os.makedirs(aligned_dir)
    
    log_align_dir = aligned_dir + '/' + 'log/'
    if not os.path.exists(log_align_dir):
        os.makedirs(log_align_dir)

    trim_read1 = trimmed1
    trim_read2 = trimmed2
    out_file = aligned_dir + '/' + sample + '_aligned.sam'

    logger_bwa_process, logger_bwa_errors = store_align_logs(log_align_dir)

    returncode = align_reads_bwa(bwa_dir, ref_fa_file, ref_index_name, trim_read1, trim_read2, 
                                                out_file, num_threads, logger_bwa_process, logger_bwa_errors)

    ##########################################################################################
    #---post_align
    ##########################################################################################
    #post_aligned_dir
    filtered_dir = out_dir + '/'+ 'filtered'
    os.makedirs(filtered_dir)
    
    log_filter_dir = filtered_dir + '/'+ 'log/'
    if not os.path.exists(log_filter_dir):
        os.makedirs(log_filter_dir)
    logger_filter_process, logger_filter_errors = store_filter_logs(log_filter_dir)
    #out_file from the align
    alignment_sam=out_file
    #min_mapq=17
    #max_soft_clip=10
    out_file1 = filtered_dir + '/' + sample + '_tmp.sam'
    stats_file = filtered_dir + '/' + sample+ '_align_stats.txt'
    primer_stats_file = filtered_dir + '/' + sample + '_primer_stats.csv'
    #max_dist = 2
    out_file2 = filtered_dir + '/' + sample + '_filtered.sam'
    filter_alignment_samtools(samtools_dir, alignment_sam, min_mapq,
                                         max_soft_clip, out_file1, stats_file,
                                         logger_filter_process, logger_filter_errors)
    identify_gs_primers(samtools_dir, out_file1, primers_file, max_dist, out_file2,
                                 primer_stats_file, stats_file, logger_filter_process,
                                 logger_filter_errors)
    ##########################################################################################
    #---Germline variant calling
    ##########################################################################################
    #germline_VC_dir
    germline_VC_dir = out_dir + '/'+ 'germline_VC'
    os.makedirs(germline_VC_dir)
    
    log_germline_VC_dir = germline_VC_dir + '/'+ 'log/'
    if not os.path.exists(log_germline_VC_dir):
        os.makedirs(log_germline_VC_dir)
    logger_germline_VC_process, logger_germline_VC_errors = store_germline_VC_logs(log_germline_VC_dir)
    
    

if __name__ == '__main__':
    main()
