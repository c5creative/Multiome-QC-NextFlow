#!/usr/bin/env nextflow

nextflow.enable.dsl=2


def get_genome (library) {
    return params.libraries[library].genome
}


def get_blacklists (genome) {
    if (params.blacklist[genome] instanceof String) {
        return [params.blacklist[genome]]
    } else {
        return params.blacklist[genome]
    }
}

process barcodes_passing_qc_thresholds {

    publishDir "${params.results}/barcodes-passing-qc-thresholds"
    tag "${library}"
    container "library://porchard/default/general:20220107"
    memory '25 GB'

    input:
    tuple val(library), path(rna), path(ataqv), path(dropkick), path(rna_barcode_list), path(atac_barcode_list)

    output:
    path("*.png")
    tuple val(library), path("${library}.pass-qc-barcodes.rna.tsv"), emit: rna_nuclei
    tuple val(library), path("${library}.pass-qc-barcodes.atac.tsv"), emit: atac_nuclei

    """
    multiome-qc.py --prefix ${library}. --dropkick $dropkick --dropkick-min ${params['libraries'][library]['thresholds']['DROPKICK']} --hqaa-min ${params['libraries'][library]['thresholds']['HQAA_MIN']} --tss-enrichment-min ${params['libraries'][library]['thresholds']['TSS_ENRICH_MIN']} --max-fraction-reads-from-single-autosome-max ${params['libraries'][library]['thresholds']['MAX_FRAC_FROM_AUTOSOME_MAX']} --umis-min ${params['libraries'][library]['thresholds']['RNA_UMI_MIN']} --fraction-mitochondrial-max ${params['libraries'][library]['thresholds']['FRAC_MITO_MAX']} $ataqv $rna $rna_barcode_list $atac_barcode_list
    """

}


process index_bam_atac {

    container 'library://porchard/default/general:20220107'

    input:
    tuple val(library), path(bam)

    output:
    tuple val(library), path(bam), path(bam_index)

    script:
    bam_index = bam.getName() + ".bai"

    """
    samtools index $bam
    """

}


process index_bam_rna {

    container 'library://porchard/default/general:20220107'

    input:
    tuple val(library), path(bam)

    output:
    tuple val(library), path(bam), path(bam_index)

    script:
    bam_index = bam.getName() + ".bai"

    """
    samtools index $bam
    """

}


// only if doing demuxlet
process clip_bam {

    memory '22 GB'
    cache 'lenient'
    tag "${library}"
    container 'library://porchard/default/general:20220107'

    input:
    tuple val(library), path(bam), path(index)

    output:
    tuple val(library), path("${library}.clipped.bam"), path("${library}.clipped.bam.bai")

    when:
    params.libraries[library].vcf != 'None'

    """
    /sw/bamUtil/bin/bam clipOverlap --poolSize 9000000 --in $bam --out ${library}.clipped.bam
    samtools index ${library}.clipped.bam
    """

}


// for RNA: do demuxlet in chunks of e.g. 500 barcodes
// ATAC is fast even w/ e.g. 10k barcodes, so don't bother for ATAC
process chunk_demuxlet_barcodes {

    cache 'lenient'
    tag "${library}"
    container 'library://porchard/default/general:20220107'

    input:
    tuple val(library), path('pass-qc-barcodes.txt')

    output:
    tuple val(library), path("${library}.barcode-batch.*.txt")

    """
    split --additional-suffix=.txt --lines 500 --suffix-length 10 pass-qc-barcodes.txt ${library}.barcode-batch.
    """

}


process demuxlet {

    container 'library://porchard/default/demuxlet:20220204'
    memory '10 GB'
    cache 'lenient'
    tag "${library} ${modality}"

    input:
    tuple val(library), path(bam), path(bam_index), val(modality), path(barcodes), path(vcf)

    output:
    tuple val(library), val(modality), path("${library}-${modality}.best"), emit: best

    """
    #demuxlet --sam $bam --vcf $vcf --alpha 0 --alpha 0.5 --group-list $barcodes --field GT --out ${library}-${modality}
    popscle demuxlet --sam $bam --vcf $vcf --alpha 0 --alpha 0.5 --group-list $barcodes --field GT --out ${library}-${modality}
    """

}


process concat_demuxlet {

    publishDir "${params.results}/demuxlet/out", mode: 'rellink'

    input:
    tuple val(library), val(modality), path("demuxlet.*.best")

    output:
    tuple val(library), val(modality), path("${library}-${modality}.best.txt")

    """
    cat demuxlet.*.best | awk 'NR==1' > header.txt
    cat header.txt > ${library}-${modality}.best.txt
    cat demuxlet.*.best | grep -v -f header.txt >> ${library}-${modality}.best.txt
    """

}


process plot_demuxlet {

    publishDir "${params.results}/demuxlet/processed", mode: 'rellink'
    container "library://porchard/default/general:20220107"

    input:
    tuple val(library), path(demuxlet_best_atac), path(demuxlet_best_rna), path(rna_barcodes), path(atac_barcodes)

    output:
    tuple val(library), path("*assignments.txt"), emit: assignments
    path("*.png")

    """
    process-demuxlet-out.py --atac-demuxlet $demuxlet_best_atac --rna-demuxlet $demuxlet_best_rna --strategy joint --prefix ${library}. --atac-barcodes $atac_barcodes --rna-barcodes $rna_barcodes
    """
}


process prep_doublet_detection {

    publishDir "${params.results}/atac-doublet-detection"
    tag "${library}"
    container "library://porchard/default/general:20220107"

    input:
    tuple val(library), path(pass_qc_barcodes), path(atac_barcode_list)

    output:
    tuple val(library), path('singlecell.csv'), path('autosomes.txt')

    """
    make-doubletdetector-files.py $pass_qc_barcodes $atac_barcode_list ${get_genome(library)}
    """

}


process run_atac_doublet_detection {

    publishDir "${params.results}/atac-doublet-detection"
    tag "${library}"
    container 'library://porchard/default/amulet:1.1'

    input:
    tuple val(library), path(bam), path(bam_index), path(single_cell), path(autosomes), path(blacklists)

    output:
    tuple val(library), path("${library}.doublet_probabilities.txt")

    """
    mkdir -p output
    zcat ${blacklists.join(' ')} | sort -k1,1 -k2n,2 > blacklist.bed
    /opt/AMULET/AMULET.sh --bcidx 0 --cellidx 1 --iscellidx 2 $bam $single_cell $autosomes blacklist.bed output/ /opt/AMULET/
    cp output/MultipletProbabilities.txt ${library}.doublet_probabilities.txt
    """

}


process plot_doublet_probabilities {

    publishDir "${params.results}/atac-doublet-detection"
    tag "${library}"
    container "library://porchard/default/general:20220107"

    input:
    tuple val(library), path(x), path(ataqv), path(rna_barcodes), path(atac_barcodes)

    output:
    path("*.png")
    path("*.singlets.txt"), emit: singlets    

    """
    plot-doublet-detector-out.py $x $ataqv $rna_barcodes $atac_barcodes $library
    """

}


workflow {

    libraries = params.libraries.keySet()

    qc_in = Channel.from(libraries.collect({it -> [it, file(params.libraries[it].rna_qc), file(params.libraries[it].atac_qc), file(params.libraries[it].rna_dropkick), file(params.rna_barcodes), file(params.atac_barcodes)]}))
    atac_bam = Channel.from(libraries.collect({it -> [it, file(params.libraries[it].atac_bam)]})) | index_bam_atac // library, bam, index
    rna_bam = Channel.from(libraries.collect({it -> [it, file(params.libraries[it].rna_bam)]})).filter({it -> params.libraries[it[0]].vcf != 'None'}) | index_bam_rna // library, bam, index

    pass_qc_metrics = barcodes_passing_qc_thresholds(qc_in)

    demuxlet_in_atac = clip_bam(atac_bam).map({it -> it + ['ATAC']}).combine(pass_qc_metrics.atac_nuclei, by: 0)
    demuxlet_rna_barcode_chunks = chunk_demuxlet_barcodes(pass_qc_metrics.rna_nuclei).transpose()
    demuxlet_in_rna = rna_bam.map({it -> it + ['RNA']}).combine(demuxlet_rna_barcode_chunks, by: 0)
    demuxlet_out = (demuxlet_in_atac.mix(demuxlet_in_rna).map({it -> it + [file(params.libraries[it[0]].vcf)]}) | demuxlet).best.groupTuple(by: [0, 1]) | concat_demuxlet // library, modality, best
    demuxlet_out_atac = demuxlet_out.filter({it -> it[1] == 'ATAC'}).map({it -> [it[0], it[2]]}) // library, best
    demuxlet_out_rna = demuxlet_out.filter({it -> it[1] == 'RNA'}).map({it -> [it[0], it[2]]}) // library, best

    dd = atac_bam.combine(prep_doublet_detection(pass_qc_metrics.atac_nuclei.combine(Channel.fromPath(params.atac_barcodes))), by: 0).map({it -> it + [get_blacklists(get_genome(it[0])).collect({x -> file(x)})]}) | run_atac_doublet_detection
    dd.map({it -> it + [file(params.libraries[it[0]].atac_qc), file(params.rna_barcodes), file(params.atac_barcodes)]}) | plot_doublet_probabilities

    demuxlet_out_atac.combine(demuxlet_out_rna, by: 0).combine(Channel.fromPath(params.rna_barcodes)).combine(Channel.fromPath(params.atac_barcodes)) | plot_demuxlet

}
