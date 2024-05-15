#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWA_MEM                 } from './modules/nf-core/bwa/mem/main'
include { IVAR_CONSENSUS          } from './modules/nf-core/ivar/consensus/main'

workflow {
    Channel
        .fromFilePairs("${params.sra_data}*_{1,2}.fastq", size: 2)
        .set { fastaq_ch }

    BWA_MEM(fastaq_ch, params.index)

    // make this a tuple of reference and gene
    Channel.from(readFastaHeaders(params.reference))
        .set { headers_ch }

    IVAR_CONSENSUS(headers_ch, BWA_MEM.out.bam)

}

def readFastaHeaders(fastaFile) {
    def headers = []
    new File(fastaFile).withReader { reader ->
        def line
        while ((line = reader.readLine()) != null) {
            if (line.startsWith(">")) {
                headers << line.substring(1)
            }
        }
    }
    return headers
}
