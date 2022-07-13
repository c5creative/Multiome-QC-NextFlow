#!/usr/bin/env python
# coding: utf-8

import sys
import glob
import json
import os
import pandas as pd

GENOME = sys.argv[1]
CELLRANGER_DIRS = sys.argv[2:]

libraries = dict()

for d in CELLRANGER_DIRS:
    library_name = os.path.basename(d)
    atac_bam = os.path.join(d, 'atac_possorted_bam.bam')
    rna_bam = os.path.join(d, 'gex_possorted_bam.bam')
    rna_matrix = os.path.join(d, 'raw_feature_bc_matrix', 'matrix.mtx.gz')
    rna_features = os.path.join(d, 'raw_feature_bc_matrix', 'features.tsv.gz')
    rna_barcodes = os.path.join(d, 'raw_feature_bc_matrix', 'barcodes.tsv.gz')
    libraries[library_name] = {
        'genome': GENOME,
        'atac_bam': atac_bam,
        'rna_bam': rna_bam,
        'rna_matrix': rna_matrix,
        'rna_features': rna_features,
        'rna_barcodes': rna_barcodes
    }

CONFIG = {'libraries': libraries}
print(json.dumps(CONFIG, sort_keys = True, indent = 4))