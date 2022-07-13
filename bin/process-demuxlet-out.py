#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from plutils import figures
from snutils import nucleus
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--atac-demuxlet', help='Demuxlet best file')
parser.add_argument('--rna-demuxlet', help='Demuxlet best file')
parser.add_argument('--atac-barcodes', help='List of whitelisted ATAC barcodes')
parser.add_argument('--rna-barcodes', help='List of whitelisted RNA barcodes')
parser.add_argument('--strategy', help='RNA, ATAC, joint (use either RNA assignment, ATAC assignment, or joint assignment)')
parser.add_argument('--prefix', default='demuxlet.', help='')
args = parser.parse_args()

PREFIX = args.prefix
RNA_BARCODES = args.rna_barcodes
ATAC_BARCODES = args.atac_barcodes

def load_demuxlet(f):
    tmp = pd.read_csv(f, sep='\t')
    return tmp

def recode_best(x):
    tmp = x.split('-')
    category = tmp[0]
    samples = '-'.join(tmp[1:])
    return f'{category} ({samples})' if category == 'SNG' else category


# read in barcodes
rna_barcodes = pd.read_csv(RNA_BARCODES, header=None)[0].to_list()
atac_barcodes = pd.read_csv(ATAC_BARCODES, header=None)[0].to_list()
if not len(rna_barcodes) == len(atac_barcodes):
    raise ValueError('Number of items in RNA barcode list != number of items in ATAC barcode list')
rna_to_atac_barcode = dict(zip(rna_barcodes, atac_barcodes))
atac_to_rna_barcode = dict(zip(atac_barcodes, rna_barcodes))


atac_demuxlet = load_demuxlet(args.atac_demuxlet).assign(modality='ATAC')
rna_demuxlet = load_demuxlet(args.rna_demuxlet).assign(modality='RNA')
atac_demuxlet['assignment'] = atac_demuxlet[['DROPLET.TYPE', 'SNG.BEST.GUESS']].apply(lambda x: '{} ({})'.format(x[0], x[1]) if x[0] == 'SNG' else x[0], axis=1)
rna_demuxlet['assignment'] = rna_demuxlet[['DROPLET.TYPE', 'SNG.BEST.GUESS']].apply(lambda x: '{} ({})'.format(x[0], x[1]) if x[0] == 'SNG' else x[0], axis=1)
assert(all(atac_demuxlet.BARCODE.isin(set(atac_barcodes))))
atac_demuxlet['RNA_BARCODE'] = atac_demuxlet.BARCODE.map(atac_to_rna_barcode)
rna_demuxlet['RNA_BARCODE'] = rna_demuxlet.BARCODE
demuxlet = atac_demuxlet[['RNA_BARCODE', 'assignment', 'NUM.SNPS']].rename(columns={'assignment': 'ATAC', 'NUM.SNPS': 'ATAC_SNPs'}).merge(rna_demuxlet[['RNA_BARCODE', 'assignment', 'NUM.SNPS']].rename(columns={'assignment': 'RNA', 'NUM.SNPS': 'RNA_SNPs'}))


demuxlet['label'] = np.where(demuxlet.RNA ==  demuxlet.ATAC, demuxlet.RNA.map(lambda x: f'both = {x}'), demuxlet[['RNA', 'ATAC']].apply(lambda x: f'RNA = {x[0]}, ATAC={x[1]}', axis=1))
# if RNA and ATAC both say singlet from same individual, assign it as such
# otherwise, assign as doublet
demuxlet['joint'] = np.where(demuxlet.RNA == demuxlet.ATAC, demuxlet.RNA, 'DBL')
demuxlet['final_assignment'] = demuxlet[args.strategy].str.split(' ').map(lambda x: x[0])
demuxlet[['RNA_BARCODE', 'final_assignment']].to_csv(f'{PREFIX}assignments.txt', sep='\t', index=False, header=False)


cmap = figures.make_colormap_dict(list(sorted(set(demuxlet.ATAC.to_list() + demuxlet.RNA.to_list()))))

summarize = demuxlet.groupby(['ATAC', 'RNA']).size().rename('n').reset_index().pivot(index='ATAC', columns='RNA', values='n').fillna(0).astype(int)
fig, ax = plt.subplots(figsize=(1+len(summarize.index), 1+len(summarize.columns)))
sns.heatmap(summarize, annot=True, ax=ax)
figures.fix_heatmap_limits(ax)
fig.savefig(f'{PREFIX}demuxlet-heatmap.png', bbox_inches='tight', dpi=300)
fig.clf()



fig, axs = plt.subplots(ncols=3, figsize=(5*3, 5))

ax = axs[0]
ax.set_title('RNA assignment')
rna_counts = demuxlet.RNA.value_counts().reset_index()
sns.barplot(x='RNA', y='index', ax=ax, data=rna_counts, palette=cmap, hue='index', dodge=False)
ax.legend().remove()
ax.set_ylabel('')
ax.set_xlabel('')

ax = axs[1]
ax.set_title('ATAC assignment')
atac_counts = demuxlet.ATAC.value_counts().reset_index()
sns.barplot(x='ATAC', y='index', ax=ax, data=atac_counts, palette=cmap, hue='index', dodge=False)
ax.legend().remove()
ax.set_ylabel('')
ax.set_xlabel('')

ax = axs[2]
ax.set_title('ATAC - RNA concordance')
counts = (demuxlet.ATAC == demuxlet.RNA).map({True: 'ATAC == RNA', False: 'ATAC != RNA'}).value_counts().rename('counts').reset_index()
sns.barplot(x='counts', y='index', ax=ax, data=counts)
ax.set_ylabel('')
ax.set_xlabel('')

fig.tight_layout()
fig.savefig(f'{PREFIX}demuxlet-bar.png', dpi=300)
fig.clf()



fig, ax = plt.subplots(figsize=(10, 8))
df = demuxlet.copy()
df.label = figures.append_n(df.label)
sns.scatterplot(x='RNA_SNPs', y='ATAC_SNPs', hue='label', ax=ax, data=df, alpha=0.1)
ax.set_xscale('log')
ax.set_yscale('log')
ax.xaxis.set_major_formatter(figures.read_count_formatter)
ax.yaxis.set_major_formatter(figures.read_count_formatter)
ax.set_ylabel('ATAC SNPs checked')
ax.set_xlabel('RNA SNPs checked')

ax.legend(bbox_to_anchor=(1.05, 1.05)).set_title('')
fig.tight_layout()
fig.savefig(f'{PREFIX}demuxlet-scatter.png', dpi=300)
fig.clf()
