#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWA_MEM                 } from '../modules/nf-core/bwa/mem/main'

workflow {
    Channel
        .fromPath('data/*.fastq.gz')
        .set { input_ch }

    BWA_MEM(input_ch)
}
