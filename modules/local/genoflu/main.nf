process GENOFLU {
    tag "$meta.id"
    label 'process_medium'
    
    publishDir path: "${params.outdir}/genoflu", mode: 'copy', enabled: false, pattern: '*.tsv'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(merged_fasta)

    output:
    tuple val(meta), path("${meta.id}.tsv"), emit: genoflu_results
    path "versions.yml",                      emit: versions

    script:
    """
    if [ ! -s ${merged_fasta} ]; then
        echo "Error: Input file ${merged_fasta} is empty or does not exist"
        exit 1
    fi

    # Capture genoflu output and format it into TSV
    genoflu.py -f ${merged_fasta} | grep "Genotype" | awk -F "--> " '{print "${meta.id}\t" \$2}' > ${meta.id}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genoflu: \$(genoflu.py --version 2>&1)
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genoflu: \$(genoflu.py --version 2>&1)
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}