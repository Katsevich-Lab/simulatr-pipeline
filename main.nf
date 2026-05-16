/*
 * simulatr-pipeline entry point.
 *
 * The actual workflow lives at workflows/simulatr.nf; main.nf just
 * imports and triggers it. This is the nf-core convention.
 */

nextflow.enable.dsl = 2

include { SIMULATR } from './workflows/simulatr.nf'

workflow {
    if (params.simulatr_specifier_fp == null) {
        log.error("--simulatr_specifier_fp is required (path to a saved simulatr_study .rds).")
        System.exit(2)
    }
    SIMULATR()
}
