# 10X Multiome QC pipeline

This pipeline subsets to quality barcodes (based on user-defined QC metric thresholds; outputs accompanying QC plot), runs demultiplexing with [demuxlet](https://github.com/statgen/popscle) (if VCF file provided), and runs ATAC-based doublet detection using [AMULET](https://github.com/UcarLab/AMULET). In addition, if starting from cellranger output or something other than our lab [snATAC-seq](https://github.com/porchard/snATACseq-NextFlow) / [snRNA-seq](https://github.com/porchard/snRNAseq-NextFlow) pipelines, this pipeline will also generate single nucleus-level RNA and ATAC QC metrics (using a custom python script and the [ataqv software](https://github.com/ParkerLab/ataqv), respectively) as well as bulk-level ATAC QC (an interactive ataqv report and, if requested, plots of bigwig signal near selected gene TSS) and bulk ATAC bigwig files, and will run [dropkick](https://github.com/KenLauLab/dropkick) on the RNA component.

Note each library must be from a single species (i.e., no barnyard experiments).

## Dependencies
[Singularity (v. 3)](https://docs.sylabs.io/guides/3.0/user-guide/) and [NextFlow](https://www.nextflow.io/) (>= v. 20.10.0). Containers with the software for each step are pulled from the Sylabs cloud library (https://cloud.sylabs.io/library) or Docker hub (https://hub.docker.com/).

## Preprocessing

If you are using the accompanying [snATAC-seq](https://github.com/porchard/snATACseq-NextFlow) and [snRNA-seq](https://github.com/porchard/snRNAseq-NextFlow) pipelines for upstream processing, you can skip this step. If you are instead starting from cellranger output or the output of some other pipeline, you'll need to run the preprocess.nf pipeline prior to running the main (qc.nf) pipeline.

The preprocessing pipeline takes the following input:
* RNA BAM file
* ATAC BAM file
* RNA UMI count matrix (market matrix format)

And produces the following:
* ATAC bigwig files
* ataqv metrics for the ATAC at the bulk and the single-nucleus levels
* RNA QC metrics
* dropkick results for the RNA

To run the pipeline, you'll need to create a `config.json` file that indicates the path to the input for each library, as well as the genome to use for each library. This should follow the format:

```python
{
    "libraries": {
        "Sample_3172-CV": {
            "atac_bam": "/lab/data/seqcore/3172-CV/10x_analysis_3172-CV/Sample_3172-CV/atac_possorted_bam.bam",
            "genome": "hg19",
            "rna_bam": "/lab/data/seqcore/3172-CV/10x_analysis_3172-CV/Sample_3172-CV/gex_possorted_bam.bam",
            "rna_barcodes": "/lab/data/seqcore/3172-CV/10x_analysis_3172-CV/Sample_3172-CV/raw_feature_bc_matrix/barcodes.tsv.gz",
            "rna_features": "/lab/data/seqcore/3172-CV/10x_analysis_3172-CV/Sample_3172-CV/raw_feature_bc_matrix/features.tsv.gz",
            "rna_matrix": "/lab/data/seqcore/3172-CV/10x_analysis_3172-CV/Sample_3172-CV/raw_feature_bc_matrix/matrix.mtx.gz"
        }
    }
}
```

You'll also need to update the file paths in the `nextflow.config` file in this directory.

Then run the pipeline:

```bin
nextflow run -resume -params-file config.json --results /path/to/results /path/to/Multiome-QC-NextFlow/preprocess.nf
```

## Running main QC pipeline

To run the main QC pipeline, you'll need to provide a config.json file like this:

```python
{
    "libraries": {
        "Sample_3172-CV-hg19": {
            "atac_bam": "/lab/work/porchard/test-multiome-qc/work/preprocess/results/atac/prune/Sample_3172-CV-hg19.pruned.bam",
            "atac_qc": "/lab/work/porchard/test-multiome-qc/work/preprocess/results/atac/ataqv/single-nucleus/Sample_3172-CV-hg19.txt",
            "genome": "hg19",
            "rna_bam": "/lab/work/porchard/test-multiome-qc/work/preprocess/results/rna/prune/Sample_3172-CV-hg19.before-dedup.bam",
            "rna_dropkick": "/lab/work/porchard/test-multiome-qc/work/preprocess/results/rna/dropkick/Sample_3172-CV-hg19.dropkick-score.tsv",
            "rna_qc": "/lab/work/porchard/test-multiome-qc/work/preprocess/results/rna/qc/Sample_3172-CV-hg19.qc.txt",
            "thresholds": { # per-nucleus QC thresholds
                "DROPKICK": "0.9", # minimum dropkick score
                "FRAC_MITO_MAX": "1.0", # maximum RNA mitochondrial contamination
                "HQAA_MIN": "6182", # minimum pass-filter ATAC read threshold
                "MAX_FRAC_FROM_AUTOSOME_MAX": "0.15", # threshold denoting the maximum fraction of filtered ATAC reads allowed to come from a single autosome; used to filter out droplets that capture e.g. a single chromosome
                "RNA_UMI_MIN": "0", # min RNA UMIs
                "TSS_ENRICH_MIN": "2" # min ATAC TSS enrichment
            },
            "vcf": "None" # path to VCF file, if demultiplexing is to be performed
        }
    }
}
```

If you've run preprocessing.nf or the accompanying ATAC/RNA pipelines, you can create this config template automatically and then manually edit any QC thresholds you wish to change:
* If you ran preprocessing.nf, you can run `python /path/to/Multiome-QC-NextFlow/bin/make-qc-config.py /path/to/preprocess/results/rna /path/to/preprocess/results/atac > config.json` to create the json
* If you ran the ATAC/RNA pipelines, you can run `python /path/to/Multiome-QC-NextFlow/bin/make-qc-config.py /path/to/rna/results /path/to/atac/results > config.json`
* Note that the `HQAA_MIN` threshold parameter will already be set to a value inferred using multi-otsu thresholding on the ATAC HQAA distribution

You'll also need to update the file paths in the `nextflow.config` file in this directory.

Then run the pipeline:

```bin
nextflow run -resume -params-file config.json --results /path/to/results /path/to/Multiome-QC-NextFlow/qc.nf
```