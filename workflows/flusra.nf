include { PROCESS_SRA             } from '../subworkflows/local/process_sra/main'
include { SRATOOLS_FASTERQDUMP    } from '../modules/nf-core/sratools/fasterqdump/main'
include { MILK_FREYJA             } from '../subworkflows/local/milk_freyja/main'

workflow FLUSRA {
    take:
    sra_accessions_ch
    milk_sra_accessions_ch
    fetch_and_pull

    main:
    sra_accessions_ch.concat(milk_sra_accessions_ch)
        .unique()
        .set { accessions_ch }

    if (params.fastq_dump_path) {
        Channel.fromFilePairs("${params.fastq_dump_path}/*_{1,2}.fastq", flat: false)
            .set { reads_ch }
    } else {
        // fetch fastq files from SRA
        SRATOOLS_FASTERQDUMP(accessions_ch)
        reads_ch = SRATOOLS_FASTERQDUMP.out.reads
    }

    if (!fetch_and_pull) {

        if (sra_accessions_ch) {
            sra_samples_ch = reads_ch
                .join(sra_accessions_ch)

            PROCESS_SRA(sra_samples_ch, params.reference)
        }

        if (milk_sra_accessions_ch) {
            milk_samples_ch = reads_ch
                .join(milk_sra_accessions_ch)

            MILK_FREYJA(milk_samples_ch, params.milk_barcode, params.milk_reference)
        }
    }
}
