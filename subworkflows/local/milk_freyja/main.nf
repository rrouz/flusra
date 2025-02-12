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
    milk_samples_ch

    main:
    MINIMAP2_ALIGN(milk_samples_ch, params.milk_reference)
    FREYJA_VARIANTS(MINIMAP2_ALIGN.out.bam, params.milk_reference)
    FREYJA_VARIANTS.out.variants.branch { meta, variants, depth ->
        to_demix: !file(depth).isEmpty()
            return tuple(meta, variants, depth)
        to_skip: true
            return tuple(meta, variants, depth)
    }.set { input_demix }
    FREYJA_DEMIX(
        input_demix.to_demix,
        params.milk_barcode,
        params.demix_autoadapt,
        params.demix_depthcutoff,
        params.demix_lineage_hierarchy
    )
}
