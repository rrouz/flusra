//
// A subworkflow providing utility functions tailored for the flusra/flusra pipeline.
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FETCH_SRA_METADATA        } from '../../../modules/local/fetch_sra_metadata'

/*
========================================================================================
    SUBWORKFLOW TO INITIALISE PIPELINE
========================================================================================
*/

workflow PIPELINE_INITIALISATION {

    take:
    bioproject
    email
    sra_metadata_file

    main:

    FETCH_SRA_METADATA(bioproject, email, sra_metadata_file)

    FETCH_SRA_METADATA.out.new_sra_metadata_file.splitText()
        .map { it.trim() }
        .set { sra_accessions_ch }

    FETCH_SRA_METADATA.out.new_milk_sra_metadata_file.splitText()
        .map { it.trim() }
        .set { milk_samples_ch }

    emit:
    sra_accessions = sra_accessions_ch
    milk_samples = milk_samples_ch
}

/*
========================================================================================
    SUBWORKFLOW FOR PIPELINE COMPLETION
========================================================================================
*/

workflow PIPELINE_COMPLETION {
    main:

    //
    // Completion summary
    //
    workflow.onComplete {
        completionSummary()
    }

    workflow.onError {
        log.error "Pipeline failed."
    }
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/
//
// Print pipeline summary on completion
//
def completionSummary() {
    if (workflow.success) {
        if (workflow.stats.ignoredCount == 0) {
            log.info "[$workflow.manifest.name]\033[0;32m Pipeline completed successfully\033[0m-"
        } else {
            log.info "[$workflow.manifest.name]\033[0;33m Pipeline completed successfully, but with errored process(es)\033[0m-"
        }
    } else {
        log.info "[$workflow.manifest.name]\033[0;31m Pipeline completed with errors\033[0m-"
    }
}
