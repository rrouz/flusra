include { BWA_MEM                 } from '../../../modules/nf-core/bwa/mem/main'
include { IVAR_CONSENSUS          } from '../../../modules/nf-core/ivar/consensus/main'
include { IVAR_VARIANTS           } from '../../../modules/nf-core/ivar/variants/main'
include { SAMTOOLS_DEPTH          } from '../../../modules/nf-core/samtools/depth/main'
include { SRATOOLS_FASTERQDUMP    } from '../../../modules/nf-core/sratools/fasterqdump/main'
include { CONSENSUS_MERGED        } from '../../../modules/local/consensus_merged/main'
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
            def parts = consensus_file.getName().split('_')
            def sampleId = parts[0]
            def gene = parts[1]
            [sampleId, gene, consensus_file]
        }
        .groupTuple(by: 0)
        .map { sampleId, genes, files -> 
            [[id: sampleId], genes, files]
        }
        .set { consensus_for_merge_ch }

    CONSENSUS_MERGED(consensus_for_merge_ch)
    
    GENOFLU(CONSENSUS_MERGED.out.merged_fasta)

    GENOFLU.out.genoflu_results
        .map { meta, tsv -> tsv }
        .collect()
        .set { genoflu_files_to_merge }

    MERGE_GENOFLU_RESULTS(genoflu_files_to_merge)

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
