include { BWA_MEM                 } from '../../../modules/nf-core/bwa/mem/main'
include { IVAR_CONSENSUS          } from '../../../modules/nf-core/ivar/consensus/main'
include { IVAR_VARIANTS           } from '../../../modules/nf-core/ivar/variants/main'
include { SAMTOOLS_DEPTH          } from '../../../modules/nf-core/samtools/depth/main'
include { GENOFLU                 } from '../../../modules/local/genoflu/main'
include { MERGE_GENOFLU_RESULTS   } from '../../../modules/local/merge_genoflu_results/main'

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
    
    IVAR_CONSENSUS.out.consensus
    .map { consensus_file -> 
        def sampleId = consensus_file.baseName.tokenize('_')[0]
        return tuple(sampleId, consensus_file)
    }
    .collectFile(
        name: { sampleId, files -> "${sampleId}.fa" },
        newLine: true
    )
    .set { consensus_for_genoflu_ch }    

    GENOFLU(consensus_for_genoflu_ch)

    GENOFLU.out.genoflu_results
    .collect()
    .set { genoflu_files_to_merge }

    // Create a channel with previous results file or a placeholder then append new results
    Channel.fromPath("${params.outdir}/genoflu/${params.genoflu_results}")
        .ifEmpty { file("NO_FILE") }
        .set { existing_results_ch }

    // Pass both new genoflu results and any existing results
    MERGE_GENOFLU_RESULTS(genoflu_files_to_merge, existing_results_ch)

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
    ).set { ch_ivar_variants_input }
    
    IVAR_VARIANTS(
        ch_ivar_variants_input,
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
