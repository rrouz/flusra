#!/usr/bin/env nextflow

include { FLUSRA  } from './workflows/flusra'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_flusra_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_flusra_pipeline'

workflow {
    if (params.bioproject) {
        PIPELINE_INITIALISATION(params.bioproject, params.email, params.metadata)
    } else {
        log.info("Skipping BioProject fetch as no BioProject ID provided")
    }

    if (!params.only_fetch) {
        if (params.bioproject && PIPELINE_INITIALISATION.out.sra_accessions) {
            FLUSRA(PIPELINE_INITIALISATION.out.sra_accessions, PIPELINE_INITIALISATION.out.milk_samples, params.fetch_and_pull)
        } else if (params.sra_accessions || params.milk_sra_accessions) {
            Channel.fromPath(params.sra_accessions)
                .splitText()
                .map { it.trim() }
                .set { sra_accessions_ch }

            Channel.fromPath(params.milk_sra_accessions)
                .splitText()
                .map { it.trim() }
                .set { milk_sra_accessions_ch }

            sra_accessions_ch.concat(milk_sra_accessions_ch)
                .unique()
                .ifEmpty {
                    log.error("No SRA accessions provided or milk samples provided")
                    exit 1
                }
                .view()

            FLUSRA(sra_accessions_ch, milk_sra_accessions_ch, params.fetch_and_pull)
        } else {
            log.info("No additional SRA accessions to process")
        }
    } else {
        log.info("Skipping SRA download and processing")
    }
    PIPELINE_COMPLETION()
}
