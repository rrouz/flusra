include { BWA_MEM                 } from '../../../modules/nf-core/bwa/mem/main'
include { IVAR_CONSENSUS          } from '../../../modules/nf-core/ivar/consensus/main'
include { IVAR_VARIANTS           } from '../../../modules/nf-core/ivar/variants/main'
include { SAMTOOLS_DEPTH          } from '../../../modules/nf-core/samtools/depth/main'


workflow PROCESS_SRA {
	take:
	sra_samples_ch

	main:
    BWA_MEM(sra_samples_ch, params.reference)

    // Generate a tuple of genes from the reference fasta file
	Channel.from(readFastaHeaders(params.reference))
            .set { genes_ch }

    IVAR_CONSENSUS(
        BWA_MEM.out.bam,
        genes_ch,
        params.reference,
        params.consensus_threshold,
        params.consensus_min_depth
    )

    SAMTOOLS_DEPTH(
        BWA_MEM.out.bam,
        genes_ch,
        params.reference
    )

    BWA_MEM.out.bam.combine(
        genes_ch.map { header ->
            def gene = header.split("\\|")[0]
            [
                gene,
                header,
                !params.gff_files.isEmpty() ? params.gff_files?.get(gene) ?: "${projectDir}/assets/NO_FILE" : "${projectDir}/assets/NO_FILE"
            ]
        }
    ).set { ivar_variants_input_ch }
    
    IVAR_VARIANTS(
        ivar_variants_input_ch,
        params.reference,
        params.variant_threshold,
        params.variant_min_depth
    )
}

def readFastaHeaders(fastaFile) {
    new File(fastaFile).readLines()
        .findAll { it.startsWith(">") }
        .collect { it.substring(1) }
    }
