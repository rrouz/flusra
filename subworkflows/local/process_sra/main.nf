include { BWA_MEM                 } from '../../../modules/nf-core/bwa/mem/main'
include { IVAR_CONSENSUS          } from '../../../modules/nf-core/ivar/consensus/main'
include { IVAR_VARIANTS           } from '../../../modules/nf-core/ivar/variants/main'
include { SAMTOOLS_DEPTH          } from '../../../modules/nf-core/samtools/depth/main'
include { SRATOOLS_FASTERQDUMP    } from '../../../modules/nf-core/sratools/fasterqdump/main'


workflow PROCESS_SRA {
	take:
	sra_samples_ch

	main:
	BWA_MEM(sra_samples_ch, params.reference)

	// Generate a tuple of params.reference and gene
	Channel.from(readFastaHeaders(params.reference))
            .set { headers_ch }

    IVAR_CONSENSUS(headers_ch, BWA_MEM.out.bam, params.reference, params.consensus_threshold, params.consensus_min_depth)
    IVAR_VARIANTS(headers_ch, BWA_MEM.out.bam, params.reference, params.gff_files, params.variant_threshold, params.variant_min_depth)
    SAMTOOLS_DEPTH(headers_ch, BWA_MEM.out.bam, params.reference)
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
