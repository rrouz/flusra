include { PROCESS_SRA             } from '../subworkflows/local/process_sra/main'
include { FASTP                   } from '../modules/nf-core/fastp/main'
include { SRATOOLS_FASTERQDUMP    } from '../modules/nf-core/sratools/fasterqdump/main'
include { SRATOOLS_PREFETCH       } from '../modules/nf-core/sratools/prefetch/main'
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
        SRATOOLS_PREFETCH(samples_ch)
        SRATOOLS_FASTERQDUMP(SRATOOLS_PREFETCH.out.sra)
        reads_ch = SRATOOLS_FASTERQDUMP.out.reads
    }

    reads_ch.branch { meta, fastq ->
        doTrim: meta.trimming_flag
            return tuple(meta, fastq,
                meta.trimming_flag.front1,
                meta.trimming_flag.front2,
                meta.trimming_flag.tail1,
                meta.trimming_flag.tail2
            )
        noTrim: true
            return tuple(meta, fastq)
    }.set { branchedFastqCh }

    branchedFastqCh.doTrim | FASTP

    output_ch = branchedFastqCh.noTrim.mix(FASTP.out.trimmed_reads)

    output_ch.multiMap { meta, reads ->
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
