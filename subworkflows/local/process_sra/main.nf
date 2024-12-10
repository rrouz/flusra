include { BWA_MEM as BWA_MEM_SRA  } from '../../../modules/nf-core/bwa/mem/main'
include { IVAR_CONSENSUS          } from '../../../modules/nf-core/ivar/consensus/main'
include { IVAR_VARIANTS           } from '../../../modules/nf-core/ivar/variants/main'
include { SAMTOOLS_DEPTH          } from '../../../modules/nf-core/samtools/depth/main'
include { SRATOOLS_FASTERQDUMP    } from '../../../modules/nf-core/sratools/fasterqdump/main'


workflow PROCESS_SRA {
	take:
	sra_samples_ch
	reference

	main:
	BWA_MEM_SRA(sra_samples_ch, reference)

	// Generate a tuple of reference and gene
	Channel.from(readFastaHeaders(reference))
            .set { headers_ch }

    IVAR_CONSENSUS(headers_ch, BWA_MEM_SRA.out.bam, reference)
    IVAR_VARIANTS(headers_ch, BWA_MEM_SRA.out.bam, reference)
    SAMTOOLS_DEPTH(headers_ch, BWA_MEM_SRA.out.bam, reference)
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
