#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// In this case, need an RNA bam, ATAC bam, and RNA matrix/barcodes/features files

IONICE = 'ionice -c2 -n7'

// Generic data
// TODO: Generic data might be best kept in the nextflow.config file
AUTOSOMAL_REFERENCES = ['hg19': (1..22).collect({it -> 'chr' + it}),
    'hg38': (1..22).collect({it -> 'chr' + it}),
    'rn4': (1..20).collect({it -> 'chr' + it}),
    'rn5': (1..20).collect({it -> 'chr' + it}),
    'rn6': (1..20).collect({it -> 'chr' + it}),
    'rn7': (1..20).collect({it -> 'chr' + it}),
    'mm9': (1..19).collect({it -> 'chr' + it}),
    'mm10': (1..19).collect({it -> 'chr' + it}),
    'mm39': (1..19).collect({it -> 'chr' + it})
]

ORGANISMS = ['hg19': 'human', 
    'hg38': 'human',
    'rn4': 'rat',
    'rn5': 'rat',
    'rn6': 'rat',
    'rn7': 'rat',
    'mm9': 'mouse',
    'mm10': 'mouse',
    'mm39': 'mouse'
]

MACS2_GENOME_SIZE = [
    'rn4': 'mm',
    'rn5': 'mm',
    'rn6': 'mm',
    'rn7': 'mm',
    'mm9': 'mm',
    'mm10': 'mm',
    'mm39': 'mm',
    'hg19': 'hs',
    'hg38': 'hs'
]


def has_blacklist (genome) {
    return params.blacklist.containsKey(genome)
}


def get_blacklists (genome) {
    if (params.blacklist[genome] instanceof String) {
        return [params.blacklist[genome]]
    } else {
        return params.blacklist[genome]
    }
}

def get_tss (genome) {
    return params.tss[genome]
}


def get_organism (genome) {
    return ORGANISMS[genome]
}


def get_chrom_sizes (genome) {
    return params.chrom_sizes[genome]
}


def get_macs2_genome_size (genome) {
    return MACS2_GENOME_SIZE[genome]
}



// RNA
process rna_prune {

    publishDir "${params.results}/rna/prune"
    maxForks 10
    tag "${library}-${genome}"
    container 'library://porchard/default/general:20220107'

    input:
    tuple val(library), val(genome), path(bam)

    output:
    tuple path("${library}-${genome}.before-dedup.bam"), path("${library}-${genome}.before-dedup.bam.bai")

    """
    ${IONICE} samtools view -h -b -q 255 -F 4 -F 256 -F 2048 $bam > ${library}-${genome}.before-dedup.bam && samtools index ${library}-${genome}.before-dedup.bam
    """

}


process rna_qc {

    memory '25 GB'
    publishDir "${params.results}/rna/qc"
    tag "${library}-${genome}"
    container 'library://porchard/default/general:20220107'

    input:
    tuple val(library), val(genome), path("star.bam"), path(matrix), path(barcodes)

    output:
    tuple val(library), val(genome), path("${library}-${genome}.qc.txt")

    """
    qc-from-starsolo.py star.bam $matrix $barcodes > ${library}-${genome}.qc.txt
    """

}


process dropkick {

    memory '150 GB'
    cpus 5
    publishDir "${params.results}/rna/dropkick"
    tag "${library}-${genome}"
    container 'library://porchard/default/dropkick:20220225'
    errorStrategy 'ignore'

    input:
    tuple val(library), val(genome), path(matrix), path(features), path(barcodes)

    output:
    tuple val(library), val(genome), path("${library}-${genome}.dropkick-score.tsv"), emit: dk_score
    path("*.png")

    """
    mkdir -p dropkick-in
    #if file $matrix | grep -q gzip; then zcat $matrix > dropkick-in/matrix.mtx; else cp $matrix dropkick-in/matrix.mtx; fi
    #if file $features | grep -q gzip; then zcat $features > dropkick-in/genes.tsv; else cp $features dropkick-in/genes.tsv; fi
    #if file $barcodes | grep -q gzip; then zcat $barcodes > dropkick-in/barcodes.tsv; else cp $barcodes dropkick-in/barcodes.tsv; fi
    if [[ $matrix == *.gz ]]; then zcat $matrix > dropkick-in/matrix.mtx; else cp $matrix dropkick-in/matrix.mtx; fi
    if [[ $features == *.gz ]]; then zcat $features > dropkick-in/genes.tsv; else cp $features dropkick-in/genes.tsv; fi
    if [[ $barcodes == *.gz ]]; then zcat $barcodes > dropkick-in/barcodes.tsv; else cp $barcodes dropkick-in/barcodes.tsv; fi
    run-dropkick.py dropkick-in/ ${library}-${genome}.
    """

}


process plot_rna_qc {

    memory '15 GB'
    publishDir "${params.results}/rna/qc"
    tag "${library}-${genome}"
    container 'library://porchard/default/dropkick:20220225'

    input:
    tuple val(library), val(genome), path(metrics), path(dk)

    output:
    tuple val(library), val(genome), path("${library}-${genome}.metrics.png"), path("${library}-${genome}.suggested-thresholds.tsv")

    """
    plot-rna-qc-metrics.py --prefix ${library}-${genome}. $metrics $dk
    """

}

// ATAC
process index_md_bam {

    time '12h'
    memory '5 GB'
    tag "${library} ${genome}"
    container 'library://porchard/default/general:20220107'

    input:
    tuple val(library), val(genome), path("${library}-${genome}.md.bam")

    output:
    tuple val(library), val(genome), path("${library}-${genome}.md.bam"), path("${library}-${genome}.md.bam.bai")

    """
    samtools index ${library}-${genome}.md.bam
    """

}


process prune {

    publishDir "${params.results}/atac/prune", mode: 'rellink', overwrite: true
    memory '3 GB'
    time '24h'
    errorStrategy 'retry'
    maxRetries 2
    tag "${library} ${genome}"
    container 'library://porchard/default/general:20220107'

    input:
    tuple val(library), val(genome), path(md_bam), path(bam_index)

    output:
    tuple val(library), val(genome), path("${library}-${genome}.pruned.bam")

    """
    ${IONICE} samtools view -h -b -f 3 -F 4 -F 8 -F 256 -F 1024 -F 2048 -q 30 $md_bam ${AUTOSOMAL_REFERENCES[genome].join(' ')} > ${library}-${genome}.unsorted.bam && samtools sort -m 2G -o ${library}-${genome}.pruned.bam -T bam-sort -O BAM ${library}-${genome}.unsorted.bam
    """

}


process bamtobed {

    time '4h'
    maxForks 10
    tag "${library} ${genome}"
    container 'library://porchard/default/general:20220107'

    input:
    tuple val(library), val(genome), path(bam)

    output:
    tuple val(library), val(genome), path("${library}-${genome}.bed")

    """
    bedtools bamtobed -i $bam > ${library}-${genome}.bed
    """

}


process macs2 {

    publishDir "${params.results}/atac/macs2", mode: 'rellink'
    time '24h'
    tag "${library} ${genome}"
    memory { 25.GB * task.attempt }
    maxRetries 2
    errorStrategy 'retry'
    container 'library://porchard/default/general:20220107'

    input:
    tuple val(library), val(genome), path(bed)

    output:
    tuple val(library), val(genome), path("${library}-${genome}_peaks.broadPeak"), emit: peaks
    tuple val(library), val(genome), path("${library}-${genome}_treat_pileup.bdg"), emit: bedgraph

    """
    macs2 callpeak -t $bed --outdir . --SPMR -f BED -n ${library}-${genome} -g ${get_macs2_genome_size(genome)} --nomodel --shift -100 --seed 762873 --extsize 200 -B --broad --keep-dup all
    """

}

process blacklist_filter_peaks {

    publishDir "${params.results}/atac/macs2", mode: 'rellink'
    time '1h'
    tag "${library} ${genome}"
    container 'library://porchard/default/general:20220107'

    input:
    tuple val(library), val(genome), path(peaks), path(blacklists)

    output:
    path("${library}-${genome}_peaks.broadPeak.noblacklist")

    when:
    has_blacklist(genome)

    """
    bedtools intersect -a $peaks -b ${blacklists.join(' ')} -v > ${library}-${genome}_peaks.broadPeak.noblacklist
    """

}


process bigwig {

    time '24h'
    publishDir "${params.results}/atac/bigwig", mode: 'rellink'
    tag "${library} ${genome}"
    container 'library://porchard/default/general:20220107'
    memory { 20.GB * task.attempt }
    maxRetries 2
    errorStrategy 'retry'

    input:
    tuple val(library), val(genome), path(bedgraph), path(chrom_sizes)

    output:
    tuple val(genome), path("${library}-${genome}.bw")

    """
    LC_COLLATE=C sort -k1,1 -k2n,2 $bedgraph > sorted.bedgraph
    bedClip sorted.bedgraph $chrom_sizes clipped.bedgraph
    bedGraphToBigWig clipped.bedgraph $chrom_sizes ${library}-${genome}.bw
    rm sorted.bedgraph clipped.bedgraph
    """

}


process plot_signal_at_tss {

    publishDir "${params.results}/atac/bigwig/plot", mode: 'rellink', overwrite: true
    errorStrategy 'retry'
    maxRetries 1
    memory { 10.GB * task.attempt }
    tag "${genome}"
    container 'library://porchard/default/general:20220107'

    input:
    tuple val(genome), path(bw), path(tss)

    output:
    path("*.png") optional true

    """
    plot-signal-at-tss.py --genes ${params.plot_signal_at_genes.join(' ')} --tss-file $tss --bigwigs ${bw.join(' ')}
    """

}


process chunk_single_nucleus_bams {

    time '48h'
    tag "${library} ${genome}"
    container 'library://porchard/default/general:20220107'

    input:
    tuple val(library), val(genome), path(md_bam), path(bam_index)

    output:
    tuple val(library), val(genome), path("${library}-${genome}.chunk*.bam")

    """
    ${IONICE} chunk-bam-by-barcode.py $md_bam ${library}-${genome}.
    """

}


process index_chunked_single_nucleus_bams {

    time '4h'
    tag "${library} ${genome} chunk_${chunk}"
    container 'library://porchard/default/general:20220107'

    input:
    tuple val(library), val(genome), val(chunk), path(bam)

    output:
    tuple val(library), val(genome), val(chunk), path(bam), path("${bam.getName() + '.bai'}")

    """
    samtools index $bam
    """

}


process ataqv_single_nucleus {

    publishDir "${params.results}/atac/ataqv/single-nucleus/json", mode: 'rellink', overwrite: true
    errorStrategy 'retry'
    maxRetries 1
    memory { 50.GB * task.attempt }
    time '10h'
    tag "${library} ${genome}"
    container 'library://porchard/default/ataqv:1.3.0'

    input:
    tuple val(library), val(genome), val(chunk), path(md_bam), path(bam_index), path(blacklists), path(tss)

    output:
    tuple val(library), val(genome), path("${library}-${genome}.chunk_${chunk}.ataqv.json.gz"), emit: json
    path("${library}-${genome}.chunk_${chunk}.ataqv.out")

    """
    export TERM=xterm-256color && ataqv --name ${library}-${genome} ${blacklists.collect({'--excluded-region-file ' + it}).join(' ')} --ignore-read-groups --nucleus-barcode-tag CB --metrics-file ${library}-${genome}.chunk_${chunk}.ataqv.json.gz --tss-file $tss ${get_organism(genome)} $md_bam > ${library}-${genome}.chunk_${chunk}.ataqv.out
    """

}


process reformat_ataqv {

    memory { 100.GB * task.attempt }
    time '10h'
    tag "${library} ${genome}"
    container 'library://porchard/default/general:20220107'

    input:
    tuple val(library), val(genome), path(json)

    output:
    tuple val(library), val(genome), path("${library}-${genome}.txt")

    """
    extractAtaqvMetric.py --files $json > ${library}-${genome}.txt
    """

}


process concat_ataqv {

    publishDir "${params.results}/atac/ataqv/single-nucleus", mode: 'rellink', overwrite: true
    time '10h'
    tag "${library} ${genome}"
    container 'library://porchard/default/general:20220107'

    input:
    tuple val(library), val(genome), path("ataqv.*.txt")

    output:
    tuple val(library), val(genome), path("${library}-${genome}.txt")

    """
    cat ataqv.*.txt | cut -f2-4 > ${library}-${genome}.txt
    """

}


process plot_qc_metrics {

    publishDir "${params.results}/atac/ataqv/single-nucleus", mode: 'rellink', overwrite: true
    time '10h'
    tag "${library} ${genome}"
    container 'library://porchard/default/dropkick:20220225'
    memory { 10.GB * task.attempt }
    maxRetries 1
    errorStrategy 'retry'

    input:
    tuple val(library), val(genome), path(metrics)

    output:
    tuple val(library), val(genome), path("*.png")
    path("*.tsv")

    """
    plot-atac-qc-metrics.py --prefix ${library}-${genome}. $metrics
    """

}


process ataqv_bulk {

    publishDir "${params.results}/atac/ataqv/bulk", mode: 'rellink', overwrite: true
    errorStrategy 'retry'
    maxRetries 1
    memory { 5.GB * task.attempt }
    time '10h'
    tag "${library} ${genome}"
    container 'library://porchard/default/ataqv:1.3.0'

    input:
    tuple val(library), val(genome), path(md_bam), path(bam_index), path(peaks), path(blacklists), path(tss)

    output:
    tuple val(genome), path("${library}-${genome}.ataqv.json.gz"), emit: json
    path("${library}-${genome}.ataqv.out")

    """
    export TERM=xterm-256color && ataqv --name ${library}-${genome} --peak-file $peaks --ignore-read-groups --metrics-file ${library}-${genome}.ataqv.json.gz --tss-file $tss ${blacklists.collect({'--excluded-region-file ' + it}).join(' ')} ${get_organism(genome)} $md_bam > ${library}-${genome}.ataqv.out
    """

}


process ataqv_bulk_viewer {

    publishDir "${params.results}/atac/ataqv/bulk", mode: 'rellink', overwrite: true
    errorStrategy 'retry'
    maxRetries 1
    memory { 10.GB * task.attempt }
    time '1h'
    tag "${genome}"
    container 'library://porchard/default/ataqv:1.3.0'

    input:
    tuple val(genome), path(json)

    output:
    path("ataqv-viewer-${genome}")

    """
    export TERM=xterm-256color && mkarv ataqv-viewer-${genome} ${json.join(' ')}
    """

}

workflow {

    libraries = params.libraries.keySet()

    // RNA
    for_rna_qc = Channel.from(libraries.collect({it -> [it, params.libraries[it].genome, file(params.libraries[it].rna_bam), file(params.libraries[it].rna_matrix), file(params.libraries[it].rna_barcodes)]}))
    for_dropkick = Channel.from(libraries.collect({it -> [it, params.libraries[it].genome, file(params.libraries[it].rna_matrix), file(params.libraries[it].rna_features), file(params.libraries[it].rna_barcodes)]}))

    Channel.from(libraries.collect({it -> [it, params.libraries[it].genome, file(params.libraries[it].rna_bam)]})) | rna_prune
    rna_qc(for_rna_qc).combine(dropkick(for_dropkick).dk_score, by: [0, 1]) | plot_rna_qc

    // ATAC
    md_bams = Channel.from(libraries.collect({it -> [it, params.libraries[it].genome, params.libraries[it].atac_bam]})) | index_md_bam
    peak_calling = md_bams | prune | bamtobed | macs2
    blacklist_filter_peaks(peak_calling.peaks.filter({it -> has_blacklist(it[1])}).map({it -> it + [get_blacklists(it[1]).collect({x -> file(x)})]}))
    bigwig(peak_calling.bedgraph.map({it -> it + [file(get_chrom_sizes(it[1]))]})).groupTuple().map({it -> it + [file(get_tss(it[0]))]}) | plot_signal_at_tss

    chunked_sn_bams = (md_bams | chunk_single_nucleus_bams).transpose().map({it -> [it[0], it[1], it[2].getName().tokenize('.')[-2].replaceAll('chunk', ''), it[2]]}) | index_chunked_single_nucleus_bams

    ((chunked_sn_bams.map({it -> it + [get_blacklists(it[1]).collect({x -> file(x)})] + [file(get_tss(it[1]))]}) | ataqv_single_nucleus).json | reformat_ataqv).groupTuple(by: [0, 1]) | concat_ataqv | plot_qc_metrics

    // ////sn_ataqv = (ataqv_single_nucleus(chunked_sn_bams)).json | reformat_ataqv).groupTuple(by: [0, 1]) | concat_ataqv | plot_qc_metrics
    // ////sn_ataqv = (((md_bams | chunk_single_nucleus_bams).transpose().map({it -> [it[0], it[1], it[2].getName().tokenize('.')[-2].replaceAll('chunk', ''), it[2]]}) | index_chunked_single_nucleus_bams | ataqv_single_nucleus).json | reformat_ataqv).groupTuple(by: [0, 1]) | concat_ataqv | plot_qc_metrics

    ataqv_bulk(md_bams.combine(peak_calling.peaks, by: [0, 1]).map({it -> it + [get_blacklists(it[1]).collect({x -> file(x)})] + [file(get_tss(it[1]))]})).json.groupTuple() | ataqv_bulk_viewer

}