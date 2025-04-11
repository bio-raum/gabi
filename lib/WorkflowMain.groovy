//
// This file holds several functions specific to the workflow/esga.nf in the nf-core/esga pipeline
//
class WorkflowMain {

    //
    // Check and validate parameters
    //
    //
    // Validate parameters and print summary to screen
    //
    public static void initialise(workflow, params, log) {
        log.info header(workflow)
    }

    public static String header(workflow) {
        def headr = ''
        def infoLine = "${workflow.manifest.description} | version ${workflow.manifest.version}"
        headr = """
    ===============================================================================
    ${infoLine}
    ===============================================================================
    """
        return headr
    }

}
