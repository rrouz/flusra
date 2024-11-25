//
// A subworkflow providing functions tailored for the analysis of milk samples using Freyja for the flusra/flusra pipeline.
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FREYJA_VARIANTS         } from '../../../modules/nf-core/freyja/variants'
include { FREYJA_DEMIX            } from '../../../modules/nf-core/freyja/demix'
include { BWA_MEM                 } from '../../../modules/nf-core/bwa/mem/main'

workflow MILK_FREYJA {
    take:
    reads
    barcode
    reference

    main:
    BWA_MEM(reads, reference)
    FREYJA_VARIANTS(BWA_MEM.out.bam, reference)
    FREYJA_DEMIX(FREYJA_VARIANTS.out.variants, barcode)
}
