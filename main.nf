#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/**
===============================
GABI Pipeline
===============================

This Pipeline performs assembly of bacterial isolates from NGS reads and performs typing

### Homepage / git
git@github.com:bio-raum/gabi.git

**/

// Pipeline version
params.version = workflow.manifest.version

include { GABI }                from './workflows/gabi'
include { BUILD_REFERENCES }    from './workflows/build_references'
include { paramsSummaryLog }    from 'plugin/nf-schema'


workflow {

    // Print summary of supplied parameters
    log.info paramsSummaryLog(workflow)


    multiqc_report = Channel.from([])
    if (!workflow.containerEngine) {
        log.warn "NEVER USE CONDA FOR PRODUCTION PURPOSES!"
    }

    WorkflowMain.initialise(workflow, params, log)
    WorkflowPipeline.initialise(params, log)

    if (params.build_references) {
        BUILD_REFERENCES()
    } else {
        GABI()
        multiqc_report = multiqc_report.mix(GABI.out.qc).toList()
    }

}
