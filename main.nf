#!/usr/bin/env nextflow

include { FLUSRA  } from './workflows/flusra'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_flusra_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_flusra_pipeline'

workflow {
    PIPELINE_INITIALISATION(params.bioproject, params.email, params.metadata)

    if (PIPELINE_INITIALISATION.out.sra_accessions) {
        FLUSRA(PIPELINE_INITIALISATION.out.sra_accessions)
    } else {
        log.info("No additional SRA accessions to process")
    }
    PIPELINE_COMPLETION()
}

