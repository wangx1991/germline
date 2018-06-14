from __future__ import barry_as_FLUFL

__all__  =  ['log_dir']
__version__  =  '1.0'
__author__  =  'Wang Xian'

import logging
import os
#import re
import sys
import time
#import gzip
#import itertools

def setup_logger(name, log_file, formatter, level=logging.DEBUG):
    handler = logging.FileHandler(log_file)
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger

def store_pipeline_logs(log_dir):
    formatter_pipeline_process = logging.Formatter("%(asctime)s;%(message)s")
    formatter_pipeline_errors = logging.Formatter(
        "%(asctime)s;%(levelname)s;                                             %(message)s")
    logger_pipeline_process = setup_logger('Running Messages of pipeline',
                                       log_dir + '/process.log',
                                       formatter_pipeline_process)
    logger_pipeline_errors = setup_logger('Errors & Warnings of pipeline',
                                      log_dir + '/errors.log',
                                      formatter_pipeline_errors)
    return logger_pipeline_process, logger_pipeline_errors

def store_trim_logs(log_dir):
    formatter_trim_process = logging.Formatter("%(asctime)s;%(message)s")
    formatter_trim_errors = logging.Formatter(
        "%(asctime)s;%(levelname)s;                                             %(message)s")
    logger_trim_process = setup_logger('Running Messages of read trimming',
                                       log_dir + '/trim_process.log',
                                       formatter_trim_process)
    logger_trim_errors = setup_logger('Errors & Warnings of read trimming',
                                      log_dir + '/trim_errors.log',
                                      formatter_trim_errors)
    return logger_trim_process, logger_trim_errors

def store_filter_logs(log_dir):
    formatter_filter_process = logging.Formatter("%(asctime)s;%(message)s")
    formatter_filter_errors = logging.Formatter("%(asctime)s;%(levelname)s; %(message)s")
    logger_filter_process = setup_logger('Post-alignment Filtration Messages',
                                         log_dir + '/filter_process.log',
                                         formatter_filter_process)
    logger_filter_errors = setup_logger('Errors & Warnings of filtration',
                                        log_dir + '/filter_errors.log',
                                        formatter_filter_errors)
    return logger_filter_process, logger_filter_errors


def store_align_logs(log_dir):
    formatter_bwa_process = logging.Formatter("%(asctime)s;%(message)s")
    formatter_bwa_errors = logging.Formatter("%(asctime)s;%(levelname)s;                                             %(message)s")
    logger_bwa_process = setup_logger('BWA Running Messages', 
                                      log_dir + '/bwa_process.log', 
                                      formatter_bwa_process)
    logger_bwa_errors = setup_logger('Errors & Warnings of BWA', 
                                    log_dir + '/bwa_errors.log',
                                    formatter_bwa_errors)
    return logger_bwa_process, logger_bwa_errors

def store_germline_VC_logs(log_dir):
    formatter_germline_VC_process = logging.Formatter("%(asctime)s;%(message)s")
    formatter_germline_VC_errors = logging.Formatter("%(asctime)s;%(levelname)s; %(message)s")
    logger_germline_VC_process = setup_logger('Germline variant calling Messages',
                                         log_dir + '/germline_VC_process.log',
                                         formatter_germline_VC_process)
    logger_germline_VC_errors = setup_logger('Errors & Warnings of Germline variant calling',
                                        log_dir + '/germline_VC_errors.log',
                                        formatter_germline_VC_errors)
    return logger_germline_VC_process, logger_germline_VC_errors