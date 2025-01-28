process FASTP {
    tag "$meta.id"
    label 'process_low', 'process_high_memory'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(reads), val(trim_front_read_1), val(trim_front_read_2), val(trim_tail_read_1), val(trim_tail_read_2)

    output:
    tuple val(meta), path('trimmed_*.fastq'), emit: trimmed_reads
    path "versions.yml"            , emit: versions

    script:
    def args = [
        trim_front_read_1 ? "--trim_front1 ${trim_front_read_1}" : null,
        trim_front_read_2 ? "--trim_front2 ${trim_front_read_2}" : null,
        trim_tail_read_1  ? "--trim_tail1 ${trim_tail_read_1}"   : null,
        trim_tail_read_2  ? "--trim_tail2 ${trim_tail_read_2}"   : null
    ].grep()
    """
    fastp \\
        -i ${reads[0]} \\
        -I ${reads[1]} \\
        ${args.join(' ')} \\
        -o trimmed_${meta.id}_1.fastq \\
        -O trimmed_${meta.id}_2.fastq

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_1.fastq
    touch ${meta.id}_2.fastq
    touch trimmed_${meta.id}_1.fastq
    touch trimmed_${meta.id}_2.fastq

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
    """
}
