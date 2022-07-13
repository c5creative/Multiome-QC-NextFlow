#!/usr/bin/env python
# coding: utf-8

import sys
import glob
import json
import os
import pandas as pd

RNA_RESULTS, ATAC_RESULTS = sys.argv[1:]

RNA_QC_FILES = glob.glob(os.path.join(RNA_RESULTS, 'qc', '*.qc.txt'))
ATAC_QC_FILES = glob.glob(os.path.join(ATAC_RESULTS, 'ataqv', 'single-nucleus', '*.txt'))
ATAC_BAMS = glob.glob(os.path.join(ATAC_RESULTS, 'prune', '*.bam'))
RNA_BAMS = glob.glob(os.path.join(RNA_RESULTS, 'prune', '*.bam'))
RNA_DROPKICK = glob.glob(os.path.join(RNA_RESULTS, 'dropkick', '*.dropkick-score.tsv'))
ATAC_SUGGESTED_HQAA = glob.glob(os.path.join(ATAC_RESULTS, 'ataqv', 'single-nucleus', '*.suggested-thresholds.tsv'))

# infer library names
rna_qc = {os.path.basename(f).replace('.qc.txt', ''): f for f in RNA_QC_FILES}
atac_qc = {os.path.basename(f).replace('.txt', ''): f for f in ATAC_QC_FILES}
atac_bam = {os.path.basename(f).replace('.pruned.bam', ''): f for f in ATAC_BAMS}
rna_bam = {os.path.basename(f).replace('.before-dedup.bam', ''): f for f in RNA_BAMS}
rna_dropkick = {os.path.basename(f).replace('.dropkick-score.tsv', ''): f for f in RNA_DROPKICK}
atac_suggested_hqaa = {os.path.basename(f).replace('.suggested-thresholds.tsv', ''): f for f in ATAC_SUGGESTED_HQAA}

# ensure that all keys are same
assert(set(rna_qc.keys()) == set(atac_qc.keys()))
assert(set(rna_qc.keys()) == set(atac_bam.keys()))
assert(set(rna_qc.keys()) == set(rna_bam.keys()))
assert(set(rna_qc.keys()) == set(rna_dropkick.keys()))
assert(set(rna_qc.keys()) == set(atac_suggested_hqaa.keys()))

# read suggested HQAA values
atac_suggested_hqaa = {l: pd.read_csv(f, sep='\t').set_index('metric').at['min_HQAA', 'threshold'] for l, f in atac_suggested_hqaa.items()}

libraries = {
    library: {
        'genome': library.split('-')[-1],
        'rna_qc': rna_qc[library],
        'atac_qc': atac_qc[library],
        'atac_bam': atac_bam[library],
        'rna_bam': rna_bam[library],
        'rna_dropkick': rna_dropkick[library],
        'vcf': 'None',
        'thresholds': {
            'HQAA_MIN': str(atac_suggested_hqaa[library]),
            'DROPKICK': '0.9',
            'TSS_ENRICH_MIN': '2',
            'MAX_FRAC_FROM_AUTOSOME_MAX': '0.15',
            'RNA_UMI_MIN': '500',
            'FRAC_MITO_MAX': '0.2'
        }
    } 
    for library in rna_qc.keys()
}



CONFIG = {'libraries': libraries}
print(json.dumps(CONFIG, sort_keys = True, indent = 4))