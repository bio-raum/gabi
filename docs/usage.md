# Usage information

Please fist check out our [installation guide](installation.md), if you haven't already done so. 

[Running the pipeline](#running-the-pipeline)

[Pipeline version](#specifying-pipeline-version)

[Choosing assembly method](#choosing-an-assembly-method)

[Options](#options)

[Expert options](#expert-options)

[Resources](#resources)

## Running the pipeline

A basic execution of the pipeline looks as follows:

a) Without a site-specific config file

```bash
nextflow run bio-raum/gabi -profile singularity --input samples.csv \\
--reference_base /path/to/references \\
--run_name pipeline-test
```

where `path_to_references` corresponds to the location in which you have [installed](installation.md) the pipeline references. 

In this example, the pipeline will assume it runs on a single computer with the singularity container engine available. Available options to provision software are:

`-profile singularity`

`-profile docker` 

`-profile podman`

`-profile conda` 

`-profile apptainer`

Additional software provisioning tools as described [here](https://www.nextflow.io/docs/latest/container.html) may also work, but have not been tested by us. Please note that conda may not work for all packages on all platforms. If this turns out to be the case for you, please consider switching to one of the supported container engines. In addition, you can set parameters such as maximum number of computing cores, RAM or the type of resource manager used (if any).

b) with a site-specific config file

```bash
nextflow run bio-raum/gabi -profile lsh --input samples.csv \\
--run_name pipeline-test 
```

In this example, both `--reference_base` and the choice of software provisioning are already set in the local configuration `lsh` and don't have to be provided as command line argument. 

## Specifying pipeline version

If you are running this pipeline in a production setting, you will want to lock the pipeline to a specific version. This is natively supported through nextflow with the `-r` argument:

```bash
nextflow run bio-raum/gabi -profile lsh -r 1.0 <other options here>
```

The `-r` option specifies a github [release tag](https://github.com/bio-raum/gabi/releases) or branch, so could also point to `main` for the very latest code release. Please note that every major release of this pipeline (1.0, 2.0 etc) comes with a new reference data set, which has the be [installed](installation.md) separately.

## Choosing an assembly method

GABI automatically chooses the appropriate assembly chain based on your data, supporting three scenarios:

- Samples with only short reads (Assembler: Shovill)
- Samples with Nanopore reads and **optional** short reads (Assembler: Dragonflye)
- Samples with only Pacbio HiFi reads (Assembler: Flye)

This is why it is important to make sure that all reads coming from the same sample are linked by a common sample ID. 

Note: HiFi data cannot be combined with any of the other technologies! (mostly because it is not necessary)

## Options

### `--input samples.csv` [default = null]

This pipeline expects a CSV-formatted sample sheet to properly pull various meta data through the processes. The required format looks as follows, depending on your input data:

#### Raw reads
If you want to assemble genomes "from scratch", you can pass raw reads:

```CSV
sample_id,platform,R1,R2
S100,ILLUMINA,/home/marc/projects/gaba/data/S100_R1.fastq.gz,/home/marc/projects/gaba/data/S100_R2.fastq.gz
```

If the pipeline sees more than one set of reads for a given sample ID and platform type, it will merge them automatically at the appropriate time. Based on what types of reads the pipeline sees, it will automatically trigger suitable tool chains. If the data set consists of only one read file (e.g. Nanopore, Pacbio), then the R2 column should remain empty. 

Please note that there is an optional column `library_id`, which is used to name some of the output folders for read-set specific QC measures. If `library_id` is not given, the pipeline will use the file name.

Allowed platforms and data types are:

* ILLUMINA (expecting PE Illumina reads in fastq format, fastq.gz)
* NANOPORE (expecting ONT reads in fastq format, fastq.gz)
* PACBIO (expecting Pacbio CCS/HiFi reads in fastq format, fastq.gz)

Read data in formats other than FastQ are not currently supported and would have to be converted into FastQ format prior to launching the pipeline. If you have a recurring use case where the input must be something other than FastQ, please let us know and we will consider it.

#### Pre-assembled genomes

You can also run GABI on pre-assembled genomes, using only those parts of the pipeline that characterize assemblies. Obviously, you will be missing many of the QC measures that rely on raw reads in one form or another. 

The required samplesheet then looks as follows:

```CSV
sample_id,assembly
S100,/path/to/S100.fasta
```

### `--build_references` [ default = null ]

This option is only used when installing the pipelines references as described [here](installation.md).

### `--fast_ref` [ default = false ]

By default, Gabi uses a comprehensive reference database to identify the best reference match per assembly. This can take a substantial amount of time, depending on completeness of the assembly and hardware. If you do not care about the best reference, but are happy with a "close enough" inference to get the correct species only, you can set this option to true. This will then run a reduced version of the database with a focus on covering relevant taxonomic groups at a much less dense sampling. Note that some of the Quast metrics may notably deteriorate as you are no longer guaranteed to get the closest possible match.

### `--run_name` [ default = null]

A name to use for various output files. This tend to be useful to relate analyses back to individual pipeline runs or projects later on. 

### `--reference_base` [ default = null ]

This option should point to the base directory in which you have installed the pipeline references. See our [installation](installation.md) instructions for details. For users who have contributed a site-specific config file, this option does not need to be set. 

### `--onthq` [ default = false ]

Set this option to true if you believe your ONT data to be of "high quality" (much of the reads >= Q20). This is typically the case for data generated with chemistry version 10.4.1 or later, preferably using a ligation protocol. This option is set to false by default.

### `--ont_min_q` [ default = 10 ]

Discard nanopore reads below this mean quality. ONT sequencing will produce a spread of qualities, typically ranging from Q10 to Q30 (the higher, the better). This option is mostly useful if you have sequenced at sufficient depth to be able to tolerate removable of some of the data in favor of higher quality reads. 

### `--ont_min_length`  [ default = 1000 ]

Discard nanopore reads below this length. Depending on your DNA extraction and/or library preparation, you will see a range of sequence lengths. If you have sequenced at sufficient depths, you may decide to discard shorter reads to improve your assembly contiguity. However, please note that discarding shorter reads may essentially throw away very short plasmids (which can be as short as ~1kb). 

## Expert options

These options are only meant for users who have a specific reason to touch them. For most use cases, the defaults should be fine. 

### `--confindr_db` [ default = null ]

A local version of the ConfindR rMLST database, available [here](https://olc-bioinformatics.github.io/ConFindr/install/#downloading-confindr-databases). Unfortunately, this database requires a personalized registration so we cannot bundle it with GABI. If no database is provided, CondindR will run without one and can consquently only use the built-in references for Escherichia, Listeria and Salmonella. 

### `--genome_size` [ default = null ]

If enabled, this is the assumed genome size against which the coverage is measured for downsampling the reads (e.g. '5Mb'). Since this pipeline supports processing of diverse species in parallel, you may wish to set this to a size that works across all expected taxa, like '6Mb'. The reads will then be downsampled to the desired max coverage, given the genome size.

### `--max_coverage` [ default = '100x']

If a genome size is specified (`--genome_size`), this is the target coverage for downsampling the read data. 

### `--max_contigs` [ default = 150 ]

If `--skip_failed` is enabled, this parameter controls the maximum number of contigs an assembly is allowed to have before it is stopped. High contig numbers are typically a sign of insufficient coverage and/or read length (in some cases it can also be a sign of excessive contamination).

### `--prokka_proteins` [ default = null ]

If you analyse a single species and wish to optimize the quality of the genome annotation, you can pass a Fasta file with known proteins to Prokka using this option, as described [here](https://github.com/tseemann/prokka?tab=readme-ov-file#option---proteins).

### `--prokka_prodigal` [ default = null ]

If you analyse a single species and wish to optimize the quality of the genome annotation, you can pass a custom prodigal training file using this option, as described [here](https://github.com/tseemann/prokka?tab=readme-ov-file#option---prodigaltf).

### `--remove_host` [ default = false ]

This option will perform filtering of short reads against a built-in reference (currently: horse) to remove any host contamination from the data. This option was found to be useful for Campylobacter, which is often grown in blood medium (in our case: horse). If you use another kind of medium and require decontamination, please open an issue and we will consider adding it. 

### `--skip_failed` [ default = false ]

By default, all samples are processed all the way to the end of the pipeline. This flag allows you to apply criteria to stop samples along the processing graph. The following criteria will be applied:

- Remove highly fragmented assemblies (see [--max_contigs](#--max_contigs))
- Remove reads that fail the ConfindR QC for intra-/inter species contamination (Illumina and Pacbio only)

### `--skip_circos` [ default = false ]

Skip generation of circos plots.

### `--skip_serotyping` [ default = false ]

Skip Serotyping

### `--skip_cgmlst` [ default = false ]

Skip cgMLST analysis

### `--skip_mlst` [ default = false ]

Skip all MLST analyses (incl. cgMLST)

### `--skip_amr` [ default = false ]

Skip prediction of AMR genes

### `--shovill_assembler` [ default = spades ]

Choose which assembly tool to use with Shovill. Valid options are skesa, velvet, megahit or spades. Default is: spades.

## Resources

The following options can be set to control resource usage outside of a site-specific [config](https://github.com/bio-raum/nf-configs) file.

### `--max_cpus` [ default = 16]

The maximum number of cpus a single job can request. This is typically the maximum number of cores available on a compute node or your local (development) machine. 

### `--max_memory` [ default = 128.GB ]

The maximum amount of memory a single job can request. This is typically the maximum amount of RAM available on a compute node or your local (development) machine, minus a few percent to prevent the machine from running out of memory while running basic background tasks.

### `--max_time`[ default = 240.h ]

The maximum allowed run/wall time a single job can request. This is mostly relevant for environments where run time is restricted, such as in a computing cluster with active resource manager or possibly some cloud environments.  
