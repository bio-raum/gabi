workflow PIPELINE_COMPLETION {

    take:
    multiqc_report

    main:
    log.info ""

}

workflow.onComplete {
    hline = '========================================='
    log.info hline
    log.info "Duration: $workflow.duration"
    log.info hline

    summary = [:]

    summary["MaxContigs"]           = params.skip_failed ? params.max_contigs : "Not applied"
    summary["ConfindR DB"]          = params.confindr_db ? params.confindr_db : "built-in"
    summary["Max Coverage"]         = params.max_coverage ? params.max_coverage : "Not applied"
    summary["Shovill assembler"]    = params.shovill_assembler
    summary["Shovill min contig length"] = params.shovill_contig_minlen
    summary["Medaka model"]         = params.medaka_model
    summary["ONT HQ reads"]         = params.onthq
    summary["ONT min Q"]            = params.ont_min_q
    summary["ONT min reads"]        = params.ont_min_reads
    summary["Homopolish model"]     = params.homopolish_model
    summary["AMRfinder"]            = [:]
    summary["Abricate"]             = [:]
    summary["AMRfinder"]["min_cov"] = params.arg_amrfinderplus_coveragemin
    summary["AMRfinder"]["min_id"]  = params.arg_amrfinderplus_identmin
    summary["Abricate"]["min_id"]   = params.arg_abricate_minid
    summary["Abricate"]["min_cov"]  = params.arg_abricate_mincov

    emailFields = [:]
    emailFields['version'] = workflow.manifest.version
    emailFields['session'] = workflow.sessionId
    emailFields['runName'] = run_name
    emailFields['success'] = workflow.success
    emailFields['dateStarted'] = workflow.start
    emailFields['dateComplete'] = workflow.complete
    emailFields['duration'] = workflow.duration
    emailFields['exitStatus'] = workflow.exitStatus
    emailFields['errorMessage'] = (workflow.errorMessage ?: 'None')
    emailFields['errorReport'] = (workflow.errorReport ?: 'None')
    emailFields['commandLine'] = workflow.commandLine
    emailFields['projectDir'] = workflow.projectDir
    emailFields['script_file'] = workflow.scriptFile
    emailFields['launchDir'] = workflow.launchDir
    emailFields['user'] = workflow.userName
    emailFields['Pipeline script hash ID'] = workflow.scriptId
    emailFields['manifest'] = workflow.manifest
    emailFields['summary'] = summary

    email_info = ''
    emailFields.each { s ->
        email_info += "\n${s.key}: ${s.value}"
    }
    params.each { p -> 
        email_info += "\n${p.key}: ${p.value}"
    }

    outputDir = new File("${params.outdir}/pipeline_info/")
    if (!outputDir.exists()) {
        outputDir.mkdirs()
    }

    outputTf = new File(outputDir, 'pipeline_report.txt')
    outputTf.withWriter { w -> w << email_info }

    // make txt template
    engine = new groovy.text.GStringTemplateEngine()

    tf = new File("$baseDir/assets/email_template.txt")
    txtTemplate = engine.createTemplate(tf).make(emailFields)
    emailText = txtTemplate.toString()

    // make email template
    hf = new File("$baseDir/assets/email_template.html")
    htmlTemplate = engine.createTemplate(hf).make(emailFields)
    emailHtml = htmlTemplate.toString()

    subject = "Pipeline finished ($run_name)."

    if (params.email) {
        mqcReport = null
        try {
            if (workflow.success && !params.skip_multiqc) {
                mqcReport = multiqc_report.getVal()
                if (mqcReport.getClass() == ArrayList) {
                    log.warn "[bio-raum/gabi] Found multiple reports from process 'multiqc', will use only one"
                    mqcReport = mqcReport[0]
                }
            }
        } catch (all) {
            log.warn '[bio-raum/gabi] Could not attach MultiQC report to summary email'
        }

        smailFields = [ email: params.email, subject: subject, emailText: emailText,
            emailHtml: emailHtml, baseDir: "$baseDir", mqcFile: mqcReport,
            mqcMaxSize: params.maxMultiqcEmailFileSize.toBytes()
        ]
        sf = new File("$baseDir/assets/sendmailTemplate.txt")
        sendmailTemplate = engine.createTemplate(sf).make(smailFields)
        sendmailHtml = sendmailTemplate.toString()

        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmailHtml
        } catch (all) {
            // Catch failures and try with plaintext
            [ 'mail', '-s', subject, params.email ].execute() << emailText
        }
    }
}