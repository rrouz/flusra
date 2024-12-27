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
include { MINIMAP2_ALIGN          } from '../../../modules/nf-core/minimap2/align'

workflow MILK_FREYJA {
    take:
    reads
    barcode
    reference

    main:
    MINIMAP2_ALIGN(reads, reference)
    FREYJA_VARIANTS(MINIMAP2_ALIGN.out.bam, reference)
    FREYJA_DEMIX(FREYJA_VARIANTS.out.variants, barcode)
}
