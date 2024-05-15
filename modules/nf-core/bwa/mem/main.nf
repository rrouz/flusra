process BWA_MEM {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta) , path(reads)
    tuple val(meta2), path(index)
    tuple val(meta3), path(fasta)

    output:
    tuple val(meta), path("*.bam")  , emit: bam,    optional: true
    path  "versions.yml"            , emit: versions

    script:
    def prefix = "${meta.id}"
    """
    bwa mem \\
        -t $task.cpus \\
        -P $index \\
        $reads \\
        | samtools view \\
        -F 4 -b \\
        | samtools sort -o ${n}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def samtools_command = sort_bam ? 'sort' : 'view'
    def extension = args2.contains("--output-fmt sam")   ? "sam" :
                    args2.contains("--output-fmt cram")  ? "cram":
                    sort_bam && args2.contains("-O cram")? "cram":
                    !sort_bam && args2.contains("-C")    ? "cram":
                    "bam"
    """
    touch ${prefix}.${extension}
    touch ${prefix}.csi
    touch ${prefix}.crai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
