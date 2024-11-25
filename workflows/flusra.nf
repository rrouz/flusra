include { BWA_MEM                 } from '../modules/nf-core/bwa/mem/main'
include { IVAR_CONSENSUS          } from '../modules/nf-core/ivar/consensus/main'
include { IVAR_VARIANTS           } from '../modules/nf-core/ivar/variants/main'
include { SAMTOOLS_DEPTH          } from '../modules/nf-core/samtools/depth/main'
include { SRATOOLS_FASTERQDUMP    } from '../modules/nf-core/sratools/fasterqdump/main'
include { MILK_FREYJA             } from '../subworkflows/local/milk_freyja/main'

workflow FLUSRA {
    take:
    sra_accessions_ch
    milk_samples_ch
    fetch_and_pull

    main:
    SRATOOLS_FASTERQDUMP(sra_accessions_ch)
    if (!fetch_and_pull) {

        BWA_MEM(SRATOOLS_FASTERQDUMP.out.reads, params.reference)

        // make this a tuple of reference and gene
        Channel.from(readFastaHeaders(params.reference))
            .set { headers_ch }

        IVAR_CONSENSUS(headers_ch, BWA_MEM.out.bam, params.reference)
        IVAR_VARIANTS(headers_ch, BWA_MEM.out.bam, params.reference)
        SAMTOOLS_DEPTH(headers_ch, BWA_MEM.out.bam, params.reference)

        if (milk_samples_ch) {
            // milk_samples_ch.view()
            // filter out the milk samples from SRA_FASTERQDUMP
            // SRATOOLS_FASTERQDUMP.out.reads
            //     .filter{ sra, _ -> sra in milk_samples_ch.toList().get() }
            //     .view()

            // print  milk_samples_ch.toList().get()

            // SRATOOLS_FASTERQDUMP.out.reads.view()

            milk_samples_reads_ch = SRATOOLS_FASTERQDUMP.out.reads
                .join(milk_samples_ch)
            
            MILK_FREYJA(milk_samples_reads_ch, params.milk_barcode, params.milk_reference)
        }
    }
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
