import os
import re
import sys
import time
import gzip
import itertools
import argparse



def script_information():
    print "\nApplication: Germline variant calling of QIAseq Targeted DNA Panel\n"
    print "====================================================================="
    print "Required environment：python \  samtools \ GATK4"

parser = argparse.ArgumentParser(usage = "\n\npython %(prog)s --source --sample_name --tailname --common_seq1 --common_seq2 --output --min_read_len --bwa_dir --ref_index_name")
parser.add_argument("--source", help = "Path to input reads in FASTA format", type = str)
parser.add_argument("--sample_name", help = "the sample name of raw reads", type = str)
parser.add_argument("--output", help = "Path of output file", type = str)
parser.add_argument("--memorySize ", help = "the cutoff of Java memory", type = str, default = '4G')
parser.add_argument("--gatk_dir", help = "the install path of GATK4", type = str)
parser.add_argument("--ref_fa_file", help = "the path of ref fasta", type = str)
parser.add_argument("--exome_target", help = "the list of exome intervals", type = str, default = NULL) 
parser.add_argument("--ERC", help = "switch to running HaplotypeCaller in GVCF mode", type = str, default = NULL)
parser.add_argument("--samtools_dir", help = "the install path of samtools", type = str)
parser.add_argument("--read_filter", help = "add a read filter that deals with some problems", type = str, default = NULL)
parser.add_argument("--reduce_logs", help = "reduce the amount of chatter in the logs: --QUIET : yes", type = str, default = 'no' , choice = ['no', 'yes'])
parser.add_argument("--create_output_variant_index", help = "turn off automatic variant index creation", type = str, default = 'false' , choice = ['false', 'true'])
parser.add_argument("-v", '-version', action = 'version', version =' %(prog)s 1.0')


#---input
    source = args.source
    sample = args.sample_name
    memorySize = "-Xmx" + args.memorySize
    gatk_dir = args.gatk_dir
    samtools_dir = args.samtools_dir
    bwa_dir = args.bwa_dir
    ref_fa_file = args.ref_fa_file
    exome_target = args.exome_target
    ERC = args.ERC
    read_filter = args.read_filter
    reduce_logs = args.reduce_logs
    