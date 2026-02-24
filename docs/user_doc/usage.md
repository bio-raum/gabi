# Usage information

Please fist check out our [installation guide](installation.md), if you haven't already done so. 

## Running the pipeline

A basic execution of the pipeline looks as follows:

=== "Built-in profile"

    ```bash
    nextflow run bio-raum/gabi \
    -r main -profile apptainer \ # (1)!
    --input samples.csv \ # (2)!
    --reference_base /path/to/references \ # (3)!
    --run_name pipeline-test
    ```

    1.  In this example, the pipeline will assume it runs on a single computer with the apptainer container engine. Available options to provision software are documented in the [installation section](installation.md).
    2.  We highly recommend pinning a release number(e.g. `-r 1.0.0`) instead of using the latest commit.
    3.  `path_to_references` corresponds to the location in which you have [installed](installation.md) the pipeline references.

=== "Site-specific profile"

    ```bash
    nextflow run bio-raum/gabi \
    -profile lsh \ # (1)!
    -r main \ # (2)!
    --input samples.csv \
    --run_name pipeline-test 
    ```

    1.  In this example, both `--reference_base` and the choice of software provisioning are already set in the  configuration `lsh` and don't have to provided as command line argument. In addition, in your site-specific configuration, you can set additional site-specific parameters, such as your local resource manager, node configuration (CPU, RAM, wall time), desired cache directory for the configured package/container software etc. It is highly recommended to [set up](installation.md) such a config file.
    2.  We highly recommend pinning a release number(e.g. `-r 1.0.0`) instead of using the latest commit.

## Removing temporary data

Nextflow stores all the process data in a folder structure inside the `work` directory. All the relevant results are subsequently copied to the designated results folder (`--outdir`). The work directory is needed to resume completed or failed pipeline runs, but should be removed once you are satisified with the analysis to save space. To do so, run:

``` bash
nextflow clean -f
```

## Specifying a pipeline version

If you are running this pipeline in a production setting, you will want to lock the pipeline to a specific version. This is natively supported through nextflow with the `-r` argument:

```bash
nextflow run bio-raum/gabi -profile myprofile -r 1.0.0 <other options here>
```

The `-r` option specifies a github [release tag](https://github.com/bio-raum/gabi/releases) or branch, so could also point to `main` for the very latest code release. Please note that every major release of this pipeline (1.0, 2.0 etc) comes with a new reference data set, which has the be [installed](installation.md) separately.

If you are feeling adventurous, you can also provide a specific commit hash, for example: `-r 575faea`.

## Updating to a new release

To use a new release, you simply have to tell Nextflow to pull the last version from the source (e.g. Github):

```bash
nextflow pull bio-raum/gabi
```

Then use the `-r` argument as explained above to run the workflow using the new release.

## Running a test

This pipeline has a built-in test to quickly check that your local setup is working correctly. To run it, do:

``` bash
nextflow run bio-raum/gabi -profile myprofile,test
```

where `myprofile` can either be a site-specific config file or one of the built-in [profiles](installation.md#software-provisioning). This test requires an active internet connection to download the test data. 

## Choosing an assembly method

GABI automatically chooses the appropriate assembly chain based on your data, supporting three scenarios:

- Samples with only short reads (Assembler: Shovill)
- Samples with Nanopore reads and **optional** short reads (Assembler: Flye + Medaka + Polypolish/Homopolish)
- Samples with only Pacbio HiFi reads (Assembler: Flye + Racon)

This is why it is important to make sure that all reads coming from the same sample are linked by a common sample ID. 

Note: HiFi data cannot be combined with any of the other technologies! (mostly because it is not necessary and usually not done)

## Command-line options

### Basic options

`--input samples.tsv` [default = null]

This pipeline expects a CSV-formatted sample sheet to properly pull various meta data through the processes. The required format looks as follows, depending on your input data.

=== "Reads"

    If you want to assemble genomes "from scratch", you can pass raw reads:

    ```TSV
    sample  platform    fq1 fq2
    S100    ILLUMINA    /home/marc/projects/gaba/data/S100_R1.fastq.gz  /home/marc/projects/gaba/data/S100_R2.fastq.gz
    ```

    !!! note Single-end sequences

        To provide Single-end sequencing data, simply leave the `fq2` field empty.

    !!! note Platform

        GABI distinguishes three sequencing platforms - ILLUMINA, NANPORE and PACBIO.

    If the pipeline sees more than one set of reads for a given sample ID (i.e. from multi-lane sequencing runs or multiple platforms), it will combine them automatically at the appropriate time.

    !!! tip Automated sample sheet generation

        If you want to automatically generate sample sheets from files in a folder, you can use a companion pipeline of GABI:

        ```bash
        nextflow run bio-raum/samplesheet -profile myprofile --input /path/to/folder --platform
        ```

        The folder should contain reads with the extension .fastq.gz or .fq.gz. The pipeline will try to guess the proper grouping rules, but please check the resulting file (results/samples.tsv) for correctness. Automatic detection of sequencing platforms will be attempted, but is not guaranteed to work perfectly for long-read data. Also note that the nested folder structure used by Nanopore demultiplexing is not currently supported.

=== "Assembly"

    You can also run GABI on pre-assembled genomes, using only those parts of the pipeline that characterize assemblies. Obviously, you will be missing many of the QC measures that rely on raw reads in one form or another. 

    The required samplesheet then looks as follows:

    ```TSV
    sample  assembly
    S100    /path/to/S100.fasta
    ```

`--reference_base` [ default = null ]

:   The location where the pipeline references are installed on your system. This will typically be pre-set in your site-specific config file and is only needed when you run without one.

    See our [installation guide](installation.md) to learn how to install the references permanently on your system.

`--run_name` [ default = null]

:   A name to use for various output files. This tends to be useful to relate analyses back to individual pipeline runs or projects later on. 

`--min_contig_len` [ default = 300 ]

:   Discard contigs shorter than this from the assembly. Very short contigs generally do not add useful information to the assembly but increase the overall size and noise. Change this value at your own discretion. 

    The default value aims to include even the shortest of (known) plasmids and some phages. 

### Illumina options

Some options specific to assembling Illumina short reads. 

`--unicycler`  [ default = false ]

:   Use [Unicycler](https://github.com/rrwick/Unicycler) over [Shovill](https://github.com/tseemann/shovill) for assembly. Shovill is an excellent tool, but hasn't been updated in many years and relies on an older version of Spades. 
    
    Unicycler in turn is slower (expect 3-4x the runtime), but potentially (!) more accurate. Note that the reference metrics used in this pipeline were derived from Shovill assemblies, so it's likely that Unicycler assemblies will (slightly) trip some of the defined thresholds (such as GC content). This doesn't necessarily mean that the Unicycler assembly is bad/worse, just different and with fewer or different kinds of 'artifacts'.  

`--shovill_assembler` [ default = spades ]

:   Choose which assembly tool to use with Shovill. Valid options are skesa, velvet, megahit or spades. Default is: spades. Incompatible with: `--unicycler`

### Nanopore options

Some options specific to assembling ONT reads. 

`--medaka_model` [ default = null ]

:   The basecalling model used for ONT reads. This option is set to null by default since more recent base callers encode this information in the sequence headers and Medaka can grab it from there. 
    
    If this is not the case for your data, you can specify the appropriate model here. 

`--homopolish_model` [ default = R10]

:   Specifies the training file to use in [Homopolish](https://github.com/ythuang0522/homopolish). Valid options are `R10` and `R9`. Most people will want to use R10. 

`--onthq` [ default = false ]

:   Set this option to true if you believe your ONT data to be of "high quality" (much of the reads >= Q20, generated with Dorado SUP basecalling). This option is set to false by default.

`--ont_min_q` [ default = 10 ]

:   Discard nanopore reads below this mean quality. ONT sequencing will produce a spread of qualities, typically ranging from Q10 to Q30 (the higher, the better). 

    This option is mostly useful if you have sequenced at sufficient depth to be able to tolerate removable of some of the data in favor of higher quality reads. 

`--ont_min_length`  [ default = 1000 ]

:   Discard nanopore reads below this length. Depending on your DNA extraction and/or library preparation, you will see a range of sequence lengths.
    
    If you have sequenced at sufficient depths, you may decide to discard shorter reads to improve your assembly contiguity. However, please note that discarding shorter reads may essentially throw away very short plasmids (which can be as short as ~1kb). 

`--skip_homopolish` [ default = false ]

:   Skip polishing using [Homopolish](https://github.com/ythuang0522/homopolish) (only the Medaka consensus assembly is used). 

    Homopolish uses homologous sequences from a database to fix potential homopolymer errors; some people may not want to include such corrections in their assembly.

`--skip_porechop` [ default = true ]

:   Skip the removal of adapters from reads using [Porechop_abi](https://github.com/bonsai-team/Porechop_ABI). 

    Porechop_abi learns potential adapter sequences directly from the read data without external knowledge. This step is skipped by default since it is a) very slow and b) because recent basecallers offer a much faster trimming option so that the reads going into GABI should typically not contain adapters anymore. 

### Expert options

These options are only meant for users who have a specific reason to touch them. For most use cases, the defaults should be fine. 

`--confindr_db` [ default = null ]

:   A local version of the ConfindR rMLST database, available [here](https://olc-bioinformatics.github.io/ConFindr/install/#downloading-confindr-databases). 

    Unfortunately, this database requires a personalized registration so we cannot bundle it with GABI. If no database is provided, CondindR will run without one and can consquently only use the built-in references for Escherichia, Listeria, Salmonella and Campylobacter. 

`--fast_ref` [ default = false ]

:   By default, GABI uses a comprehensive reference database to identify the best reference match per assembly. This can take a substantial amount of time, depending on completeness of the assembly and hardware. 

    If you do not care about the best reference, but are happy with a "close enough" inference to get the correct species only, you can set this option to true. This will then run a reduced version of the database with a focus on covering relevant taxonomic groups at a much less dense sampling. Note that some of the Quast metrics may notably deteriorate as you are no longer guaranteed to get the closest possible match. This approach may yield subpar results if your sample belongs to a group of closely related taxa, such as <i>Campylobacter</i>.

`--max_coverage` [ default = null ]

:   Performs downsampling of the read data to the specified depth. This is done for each sequencing platform, so if you have both Illumina and ONT reads for a given sample, each set will be downsampled separately. 
    
    Please not that downsampling uses a random seed to choose which reads to retain and will thus yield slightly differing results each time and/or on different systems. Use with `--random_seed` if you would like the results to be reproducible.

`--max_contigs` [ default = 150 ]

:   If `--skip_failed` is enabled, this parameter controls the maximum number of contigs an assembly is allowed to have before it is stopped. 

    High contig numbers are typically a sign of insufficient coverage and/or read length (in some cases it can also be a sign of excessive contamination).

`--prokka_proteins` [ default = null ]

:   If you analyse a single species and wish to optimize the quality of the genome annotation, you can pass a Fasta file with known proteins to Prokka using this option, as described [here](https://github.com/tseemann/prokka?tab=readme-ov-file#option---proteins).

`--prokka_prodigal` [ default = null ]

:   If you analyse a single species and wish to optimize the quality of the genome annotation, you can pass a custom prodigal training file using this option, as described [here](https://github.com/tseemann/prokka?tab=readme-ov-file#option---prodigaltf).

`--random_seed` [ default = false ]

:   A random seed to use during read downsampling (when using --max_coverage). Downsampling will randomly choose reads to retain; if you need your analysis to be perfectly reproducible, provide a random number as seed to fix the sampling to a specific set of reads. 

`--remove_host` [ default = false ]

:   This option will perform filtering of short reads against a built-in reference (currently: horse) to remove any host contamination from the data. 

    This option was found to be useful for Campylobacter, which is often grown in blood medium (in our case: horse). If you use another kind of medium and require decontamination, please open an issue and we will consider adding it. 

`--skip_failed` [ default = false ]

:   By default, all samples are processed all the way to the end of the pipeline. This flag allows you to apply criteria to stop samples along the processing graph. The following criteria will be applied:

    - Remove highly fragmented assemblies
    - Remove reads that fail the ConfindR QC for intra-/inter species contamination (Illumina and Pacbio only)

`--skip_amr` [ default = false ]

:   Skip prediction of AMR genes

`--skip_annotation` [ default = false ]

:   Skip annotation of gene models

`--skip_circos` [ default = false ]

:   Skip generation of circos plots.

`--skip_mlst` [ default = false ]

:   Skip MLST analyses

`--skip_optional` [ default = false ]

:   Short-cut to skip mlst, amr, variant analyses, serotyping and annotation (equivalent to: `--skip_amr --skip_mlst --skip_variants --skip_annotation --skip_serotyping`)

`--skip_serotyping` [ default = false ]

:   Skip Serotyping

`--skip_variants` [ default = false ]

:   Skip variant calling
