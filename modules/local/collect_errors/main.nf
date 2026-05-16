process COLLECT_ERRORS {
    label 'COLLECT_ERRORS'

    input:
    path error_jsons

    publishDir params.result_dir, mode: 'copy'

    output:
    path 'errors.parquet', optional: true
    path 'errors.tsv',     optional: true

    script:
    """
    collect_errors.R errors.tsv errors.parquet
    """
}
