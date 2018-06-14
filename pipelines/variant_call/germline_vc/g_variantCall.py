from __future__ import barry_as_FLUFL

__all__  =  ['--source' , 'sample_name ' , 'output ' , 'memorySize' , 'gatk_dir ' , 'ref_fa_file' , 'exome_target' , 'ERC' , 'samtools_dir' , 'read_filter' , 'reduce_logs' , 'create_output_variant_index']
__version__  =  '1.0'
__author__  =  'Wang Xian'

import os
import re
import sys
import time
import gzip
import itertools
import argparse


