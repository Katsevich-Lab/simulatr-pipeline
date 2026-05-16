process COLLECT_MEMORY {
    label 'COLLECT_MEMORY'

    input:
    path memory_tsvs

    publishDir params.result_dir, mode: 'copy'

    output:
    path 'memory.parquet', optional: true
    path 'memory.tsv'

    script:
    """
    collect_memory.R memory.tsv memory.parquet
    """
}
