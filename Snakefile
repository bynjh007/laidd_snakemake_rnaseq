from os import path

import numpy as np
import pandas as pd

################################################################################
# Globals                                                                      #
################################################################################

samples = pd.read_csv(config['samples'])

# default path for analysis
data_path = config["path"]["default"]

# directory for raw files
raw_dir = path.join(data_path, 'raw')

# directory for bam files
bam_dir = path.join(data_path, 'bam')

# directory for samtools stat files
stats_dir = path.join(data_path, 'stats')

featureCounts_dir=path.join(data_path, 'featureCounts')

qc_dir=path.join(data_path, 'qc')
log_dir=path.join(data_path, 'log')


################################################################################
# Functions                                                                    #
################################################################################

all_samples = samples.SAMPLE_ID.tolist()

def format_options(options):
    return ' '.join(options or [])

################################################################################
# Rules                                                                        #
################################################################################

def all_input(wildcards):
    quant_out = path.join(featureCounts_dir, 'merged.gene.txt')
    qc_out = [path.join(qc_dir, 'multiqc_report.html')]
    all_list = quant_out + qc_out    
    return all_list


rule all:
    input:
        all_input

include: 'rules/cutadapt.smk'
include: 'rules/qc.smk'
include: 'rules/alignment.smk'
include: 'rules/quantification.smk'

