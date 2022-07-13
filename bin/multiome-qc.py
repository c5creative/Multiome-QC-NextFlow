#!/usr/bin/env python
# coding: utf-8

import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import re
import snutils.nucleus as nucleus
import argparse

parser = argparse.ArgumentParser(add_help=True)
parser.add_argument('atac_metrics', help='Path to file of ataqv metrics (barcode, metric, value)')
parser.add_argument('rna_metrics', help='Path to file of RNA metrics')
parser.add_argument('rna_barcodes', help='Path to file of RNA whitelisted barcodes. The order must match the order in the ATAC whitelisted barcodes file (e.g., the first ATAC barcode corresponds to the first RNA barcode, etc)')
parser.add_argument('atac_barcodes', help='Path to file of ATAC whitelisted barcodes. The order must match the order in the RNA whitelisted barcodes file (e.g., the first ATAC barcode corresponds to the first RNA barcode, etc)')
parser.add_argument('--dropkick', default=None, help='Path to file of dropkick scores (barcode, score; headered)')
parser.add_argument('--prefix', default='qc.', help='Prefix for output files (default: "qc.")')
parser.add_argument('--dropkick-min', dest='dropkick_min', type=float, default=0.0, help='Min. dropkick score (default: 0)')
parser.add_argument('--hqaa-min', dest='hqaa_min', type=int, default=20000, help='Min. HQAA to pass QC (default: 20000)')
parser.add_argument('--tss-enrichment-min', dest='tss_enrichment_min', type=float, default=2, help='Min. TSS enrichment to pass QC (default: 2)')
parser.add_argument('--max-fraction-reads-from-single-autosome-max', dest='max_frac_from_autosome_max', type=float, default=0.15, help='Max. max fraction reads from single autosome to pass QC (default: 0.15)')
parser.add_argument('--umis-min', dest='umis_min', type=int, default=500, help='Min. UMIs to pass QC (default: 500)')
parser.add_argument('--fraction-mitochondrial-max', dest='rna_frac_mito_max', type=float, default=0.02, help='Max. fraction mitochondrial UMIs to pass QC (default: 0.02)')
args = parser.parse_args()

ATAQV_METRIC_FILE = args.atac_metrics
RNAQC_METRIC_FILE = args.rna_metrics
ATAC_BARCODES = args.atac_barcodes
RNA_BARCODES = args.rna_barcodes
DROPKICK = args.dropkick
DROPKICK_MIN = args.dropkick_min
HQAA_MIN = args.hqaa_min
UMIS_MIN = args.umis_min
TSS_ENRICHMENT_MIN = args.tss_enrichment_min
MAX_FRAC_FROM_AUTOSOME_MAX = args.max_frac_from_autosome_max
RNA_FRAC_MITO_MAX = args.rna_frac_mito_max
PREFIX = args.prefix

@ticker.FuncFormatter
def read_count_formatter(x, pos):
    if x >= 1e9:
        return '{}B'.format(x/1e9)
    if x >= 1e6:
        return '{}M'.format(x/1e6)
    if x >= 1e3:
        return '{}k'.format(x/1e3)
    else:
        return x


def load_rnaqc_metric_file(f):
    tmp = pd.read_csv(f, sep='\t')
    tmp = tmp[~tmp.barcode.isin(['-', 'no_barcode'])].set_index('barcode')
    for col in tmp.columns:
        if 'alignment' in col or 'reads' in col or col == 'umis':
            tmp[col] = tmp[col].astype(int)
        elif 'fraction' in col:
            tmp[col] = tmp[col].astype(float)
    tmp = tmp.reset_index()
    return tmp


# read in barcodes
rna_barcodes = pd.read_csv(RNA_BARCODES, header=None)[0].to_list()
atac_barcodes = pd.read_csv(ATAC_BARCODES, header=None)[0].to_list()
if not len(rna_barcodes) == len(atac_barcodes):
    raise ValueError('Number of items in RNA barcode list != number of items in ATAC barcode list')
rna_to_atac_barcode = dict(zip(rna_barcodes, atac_barcodes))
atac_to_rna_barcode = dict(zip(atac_barcodes, rna_barcodes))

# read in QC metrics
ataqv_metrics = pd.read_csv(ATAQV_METRIC_FILE, sep='\t', header=None, names=['barcode', 'metric', 'value'])
ataqv_metrics = ataqv_metrics[ataqv_metrics.barcode!='no_barcode']
rnaqc_metrics = load_rnaqc_metric_file(RNAQC_METRIC_FILE)

dropkick_dict = None
if DROPKICK:
    dropkick = pd.read_csv(DROPKICK, sep='\t')
    dropkick.columns = ['barcode', 'dropkick_score']
    dropkick_dict = dict(zip(dropkick.barcode, dropkick.dropkick_score))

# combine metrics from different libraries
# use RNA barcodes
assert(all(ataqv_metrics.barcode.isin(set(atac_barcodes))))
ataqv_metrics.barcode = ataqv_metrics.barcode.map(atac_to_rna_barcode)
ataqv_metrics = ataqv_metrics.pivot(index='barcode', columns='metric', values='value').loc[:,['hqaa', 'tss_enrichment', 'max_fraction_reads_from_single_autosome']].reset_index()
rnaqc_metrics = rnaqc_metrics[['barcode', 'umis', 'fraction_mitochondrial']]

multiome_metrics = ataqv_metrics.merge(rnaqc_metrics, on=['barcode'])
multiome_metrics.hqaa = multiome_metrics.hqaa.astype(int)
multiome_metrics.umis = multiome_metrics.umis.astype(int)
multiome_metrics.tss_enrichment = multiome_metrics.tss_enrichment.map(lambda x: '0' if x == 'None' else x).astype(float)
multiome_metrics.max_fraction_reads_from_single_autosome = multiome_metrics.max_fraction_reads_from_single_autosome.astype(float)
multiome_metrics.fraction_mitochondrial = multiome_metrics.fraction_mitochondrial.astype(float)
multiome_metrics['dropkick_score'] = multiome_metrics.barcode.map(dropkick_dict) if dropkick_dict else 0

multiome_metrics['pass_qc'] = (multiome_metrics.hqaa >= HQAA_MIN) & (multiome_metrics.umis >= UMIS_MIN) & (multiome_metrics.tss_enrichment >= TSS_ENRICHMENT_MIN) & (multiome_metrics.max_fraction_reads_from_single_autosome <= MAX_FRAC_FROM_AUTOSOME_MAX) & (multiome_metrics.fraction_mitochondrial <= RNA_FRAC_MITO_MAX) & (multiome_metrics.dropkick_score >= DROPKICK_MIN)
pass_qc_barcodes_rna = multiome_metrics[multiome_metrics.pass_qc].barcode.to_list()
pass_qc_barcodes_atac = [rna_to_atac_barcode[i] for i in pass_qc_barcodes_rna]

with open(f'{PREFIX}pass-qc-barcodes.rna.tsv', 'w') as fh:
    fh.write('\n'.join(pass_qc_barcodes_rna) + '\n')

with open(f'{PREFIX}pass-qc-barcodes.atac.tsv', 'w') as fh:
    fh.write('\n'.join(pass_qc_barcodes_atac) + '\n')

# plot HQAA vs UMIs
# plot HQAA vs TSS enrichment
# plot HQAA vs max fraction autosomal
# plot UMIs vs fraction mitochondrial
fig, axs = plt.subplots(ncols=2, nrows=2, figsize=(10, 10))

number_pass_qc = sum(multiome_metrics.pass_qc==True)
number_fail_qc = sum(multiome_metrics.pass_qc==False)
recode = {True: 'True (n={:,})'.format(number_pass_qc), False: 'False (n={:,})'.format(number_fail_qc)}
label_colors = {recode[True]: 'red', recode[False]: 'black'}
multiome_metrics['label'] = multiome_metrics.pass_qc.map(recode)

# HQAA vs UMIs
ax = axs[0,0]
g = sns.scatterplot(x='umis', y='hqaa', hue='label', palette=label_colors, data=multiome_metrics, alpha=0.05, ax=ax)
g.legend().set_title('Pass QC')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(left=1)
ax.set_ylim(bottom=1)
ax.xaxis.set_major_formatter(read_count_formatter)
ax.yaxis.set_major_formatter(read_count_formatter)
ax.set_xlabel('RNA UMIs')
ax.set_ylabel('ATAC pass filter reads')
ax.axvline(UMIS_MIN, color='red', linestyle='dashed')
ax.axhline(HQAA_MIN, color='red', linestyle='dashed')

# HQAA vs TSS enrichment
ax = axs[0,1]
g = sns.scatterplot(x='hqaa', y='tss_enrichment', hue='label', palette=label_colors, data=multiome_metrics, alpha=0.05, ax=ax)
g.legend().set_title('Pass QC')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(left=10)
ax.set_ylim(top=100, bottom=0.1)
ax.xaxis.set_major_formatter(read_count_formatter)
ax.set_xlabel('ATAC pass filter reads')
ax.set_ylabel('ATAC TSS enrichment')
ax.axvline(HQAA_MIN, color='red', linestyle='dashed')
ax.axhline(TSS_ENRICHMENT_MIN, color='red', linestyle='dashed')

# HQAA vs max fraction autosomal, human
ax = axs[1,0]
g = sns.scatterplot(x='hqaa', y='max_fraction_reads_from_single_autosome', hue='label', palette=label_colors, data=multiome_metrics, alpha=0.01, ax=ax)
g.legend().set_title('Pass QC')
ax.set_xscale('log')
ax.set_xlim(left=100)
ax.xaxis.set_major_formatter(read_count_formatter)
ax.set_xlabel('ATAC pass filter reads')
ax.set_ylabel('Max. fraction reads from single autosome')
ax.axvline(HQAA_MIN, color='red', linestyle='dashed')
ax.axhline(MAX_FRAC_FROM_AUTOSOME_MAX, color='red', linestyle='dashed')

# UMIs vs fraction mitochondrial, human
ax = axs[1,1]
g = sns.scatterplot(x='umis', y='fraction_mitochondrial', hue='label', palette=label_colors, data=multiome_metrics, alpha=0.01, ax=ax)
g.legend().set_title('Pass QC')
ax.set_xscale('log')
ax.set_xlim(left=10)
ax.xaxis.set_major_formatter(read_count_formatter)
ax.set_xlabel('RNA UMIs')
ax.set_ylabel('RNA frac. mitochondrial')
ax.axvline(UMIS_MIN, color='red', linestyle='dashed')
ax.axhline(RNA_FRAC_MITO_MAX, color='red', linestyle='dashed')

fig.tight_layout()
fig.savefig(f'{PREFIX}plots.png')
fig.clf()
