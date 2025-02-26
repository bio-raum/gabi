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

        // Print help to screen if required
        if (params.help) {
            log.info help(workflow)
            System.exit(0)
        }
    }

    // TODO: Change name of the pipeline below
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

    public static String help(workflow) {
        def command = "nextflow run ${workflow.manifest.name} -r <VERSION> -profile <PROFILE> --input some_file.csv --run_name some_name"
        def helpString = ''
        // Help message
        helpString = """

            Usage: $command

            Required parameters:
            --input                        The primary pipeline input (typically a CSV file)
            --email                        Email address to send reports to (enclosed in '')
            --run_name                     A descriptive name for this pipeline run
            --reference_base               Location of the locally installed references
            Optional parameters:
             Illumina:
             --shovill_assembler           Which assembler to use by Shovill (skesa, velvet, megahit oder spades [default])
             --shovill_contig_minlen       Minimum length of contigs to keep after assembly (default: 600)
             Nanopore:
             --medaka_model                Basecalling model used to call ONT reads (default: null, autodetect)
             --homopolish_model            Model to use for Homopolish (R10 [default], R9)
             --skip_homopolish             Skip polishing of ONT assemblies with Homopolish
             --onthq                       Nanopore reads are "high quality" (v5.0.0 SUP)
             --ont_min_q                   Minimum phred score (Q) of reads to keep (default: 10)
             --ont_min_length              Minimum length of ONT reads to keep (default 1000)
            Expert options:
            --confindr_db                  Provide your own database to ConfindR (default: null)
            --fast_ref                     Identify species from assembly with a reduced database (faster, but slightly less accurate)
            --max_coverage                 When set, downsample reads to this approx. coverage (default: '100x')
            --skip_amr                     Skip prediction of AMR genes
            --skip_mlst                    Skip MLST typing
            --skip_serotyping              Skip serotyping
            --max_cpus                     Maximum amount of CPUs to use (default: 8)
            --max_memory                   Maximum amount of RAM to use (default: '64.GB')
            --max_time                     Maximum walltime to allow (HPC clusters only, default = '240.h')
            Output:
            --outdir                       Local directory to which all output is written (default: results)
        """
        return helpString
    }

}
