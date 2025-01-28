include { PROCESS_SRA             } from '../subworkflows/local/process_sra/main'
include { SRATOOLS_FASTERQDUMP    } from '../modules/nf-core/sratools/fasterqdump/main'
include { MILK_FREYJA             } from '../subworkflows/local/milk_freyja/main'

workflow FLUSRA {
    take:
    samples_ch

    main:
    if (params.fastq_dump_path) {
        // Allow for the user to provide the fastq files directly instead of fetching from SRA
        Channel.fromFilePairs("${params.fastq_dump_path}/*_{1,2}.fastq", flat: false)
            .map { id, pair -> tuple(pair, id) }
            .set { fastq_ch }

        // Join the fastq files with the samples channel
        samples_ch.join(
            fastq_ch,
            by: 1,
            failOnDuplicate: true
        ).map { sra, meta, reads -> 
            tuple(meta, reads)
        }.set { reads_ch }
    } else {
        // Fetch fastq files from SRA
        SRATOOLS_FASTERQDUMP(samples_ch)
        reads_ch = SRATOOLS_FASTERQDUMP.out.reads
    }

    reads_ch.multiMap { meta, reads ->
        samples: meta.process_flag ? tuple(meta, reads) : null
        milk: meta.milk_flag ? tuple(meta, reads) : null
    }.set { sample_reads_input }

    if (!params.fetch_and_pull) {
        sample_reads_input
            .samples
            .filter { it != null } | PROCESS_SRA

        sample_reads_input
            .milk
            .filter { it != null } | MILK_FREYJA
    }
}
