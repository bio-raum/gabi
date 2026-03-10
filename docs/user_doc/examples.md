# Examples running GABI

Below are *some* examples of how you could approach assembly of your genomes. Below, we will assume that you use a site-specific config file named `my_profile`. See our [installation](../user_doc/installation.md) for ways to configure the pipeline. 

## Illumina short-reads

This is probably the most common approach used for obtaining genomic information from bacterial isolates. Typically, you will have one set of paired-end reads per sample - so your samplesheet would look as follows:

```TSV
sample  platform    fq1 fq2
SampleA ILLUMINA    /path/to/sampleA_R1.fastq.gz    /path/to/sampleA_R2.fastq.gz
SampleB ILLUMINA    /path/to/sampleB_R1.fastq.gz    /path/to/sampleB_R2.fastq.gz
```

With this samplesheet, GABI can then be run like so:

=== "All default options"

    Run GABI with all-defaults options. 

    ```BASH
    nextflow run bio-raum/gabi -profile my_profile \
    -r 1.4.0 \
    --input samples.tsv \
    --run_name MyIlluminaRun
    ```

=== "Using Unicycler"

    Switch from Shovill to Unicycler for assembly. 

    ```BASH
    nextflow run bio-raum/gabi -profile my_profile \
    -r 1.4.0 \
    --input samples.tsv \
    --unicycler \
    --run_name MyIlluminaRun 
    ```

=== "Assembly-only"

    Only assemble genomes and do not perform downstream characterization.

    ```BASH
    nextflow run bio-raum/gabi -profile my_profile \
    -r 1.4.0 \
    --input samples.tsv \
    --skip_optional \
    --run_name MyIlluminaRun 
    ```

## ONT reads

Inexpensive long-reads generated using Nanopore sequencing are an attractive alternative to Illumina sequencing as the added information per read typically enables the reconstruction of the entire bacterial chromosome into one contig. Downsides include platform-specific homopolymer artifacts and a lower per-base quality, the latter of which can be largely compensated by added sequencing depth. 

Nanopore data per sample is typically split across many read files - which you could either merge before running GABI or provide individually with the same sample ID. With Nanopore being a single-end technology, the `fq2` remains empty, of course. 

```TSV
sample  platform    fq1 fq2
SampleA NANOPORE    /path/to/BBJ413_pass_barcode01_34ed0b48_e8967f4e_0.fastq.gz
SampleA NANOPORE    /path/to/BBJ413_pass_barcode01_34ed0b48_e8967f4e_1.fastq.gz
SampleA NANOPORE    /path/to/BBJ413_pass_barcode01_34ed0b48_e8967f4e_2.fastq.gz
SampleA NANOPORE    /path/to/BBJ413_pass_barcode01_34ed0b48_e8967f4e_3.fastq.gz
SampleB NANOPORE    /path/to/BBJ413_pass_barcode01_34ed0b48_e8967f4e_0.fastq.gz
SampleB NANOPORE    /path/to/BBJ413_pass_barcode01_34ed0b48_e8967f4e_1.fastq.gz
SampleB NANOPORE    /path/to/BBJ413_pass_barcode01_34ed0b48_e8967f4e_2.fastq.gz
SampleB NANOPORE    /path/to/BBJ413_pass_barcode01_34ed0b48_e8967f4e_3.fastq.gz
```

With this samplesheet, you can run GABI:

=== "All default options"

    Run GABI with all-defaults options. 

    ```BASH
    nextflow run bio-raum/gabi -profile my_profile \
    -r 1.4.0 \
    --input samples.tsv \
    --run_name MyIlluminaRun
    ```

=== "Consensus assembly"

    Perform consensus assembly with Autocycler.

    ```BASH
    nextflow run bio-raum/gabi -profile my_profile \
    -r 1.4.0 \
    --input samples.tsv \
    --autocycler \
    --run_name MyIlluminaRun
    ```

=== "Consensus assembly with homopolish"

    Perform consensus assembly with Autocycler and polish with Homopolish.

    ```BASH
    nextflow run bio-raum/gabi -profile my_profile \
    -r 1.4.0 \
    --input samples.tsv \
    --autocycler \
    --homopolish \
    --run_name MyIlluminaRun
    ```

## ONT and Illumina reads

Combining Nanopore long reads with Illumina short reads will extend the assembly process by using short reads for assembly polishing. No additional command-line flags are required. A samplesheet with mixed ONT and Illumina reads will look as follows:

```TSV
sample  platform    fq1 fq2
SampleA ILLUMINA    /path/to/sampleA_R1.fastq.gz    /path/to/sampleA_R2.fastq.gz
SampleA NANOPORE    /path/to/BBJ413_pass_barcode01_34ed0b48_e8967f4e_0.fastq.gz
SampleA NANOPORE    /path/to/BBJ413_pass_barcode01_34ed0b48_e8967f4e_1.fastq.gz
SampleA NANOPORE    /path/to/BBJ413_pass_barcode01_34ed0b48_e8967f4e_2.fastq.gz
SampleA NANOPORE    /path/to/BBJ413_pass_barcode01_34ed0b48_e8967f4e_3.fastq.gz
SampleB ILLUMINA    /path/to/sampleB_R1.fastq.gz    /path/to/sampleB_R2.fastq.gz
SampleB NANOPORE    /path/to/BBJ413_pass_barcode01_34ed0b48_e8967f4e_0.fastq.gz
SampleB NANOPORE    /path/to/BBJ413_pass_barcode01_34ed0b48_e8967f4e_1.fastq.gz
SampleB NANOPORE    /path/to/BBJ413_pass_barcode01_34ed0b48_e8967f4e_2.fastq.gz
SampleB NANOPORE    /path/to/BBJ413_pass_barcode01_34ed0b48_e8967f4e_3.fastq.gz
```

Note that Homopolish will not run if short reads are provided. 

## Pacbio reads

Assembling genomes form Pacbio reads is largely equivalent to using ONT reads, with the exception that Pacbio distinguishes between corrected (HiFi) and uncorrected (CLR) reads - which GABI needs to be told about. GABI only accepts one kind (HiFI or CLR) of Pacbio reads per run. 

```TSV
sample  platform    fq1 fq2
sampleA PACBIO  /path/to/sampleA_hifi.fastq.gz
sampleB PACBIO  /path/to/sampleB_hifi.fastq.gz
```

GABI can then be run like so:

=== "CLR reads, default options"

    Run GABI with all-defaults options for CLR reads. 

    ```BASH
    nextflow run bio-raum/gabi -profile my_profile \
    -r 1.4.0 \
    --input samples.tsv \
    --run_name MyIlluminaRun
    ```

=== "HiFi reads, default options"

    Run GABI with all-defaults options for HiFi reads. 

    ```BASH
    nextflow run bio-raum/gabi -profile my_profile \
    -r 1.4.0 \
    --input samples.tsv \
    --pacbio_hifi \
    --run_name MyIlluminaRun
    ```

=== "CLR reads with Homopolish"

    Run GABI with CLR reads and use Homopolish

    ```BASH
    nextflow run bio-raum/gabi -profile my_profile \
    -r 1.4.0 \
    --input samples.tsv \
    --homopolish \
    --run_name MyIlluminaRun
    ```

## Pacbio and Illumina reads

Combining Pacbio long reads with Illumina short reads will extend the assembly process by using short reads for assembly polishing. No additional command-line flags are required. A samplesheet with mixed Pacbio and Illumina reads will look as follows:

```TSV
sample  platform    fq1 fq2
SampleA ILLUMINA    /path/to/sampleA_R1.fastq.gz    /path/to/sampleA_R2.fastq.gz
sampleA PACBIO  /path/to/sampleA_hifi.fastq.gz
SampleB ILLUMINA    /path/to/sampleB_R1.fastq.gz    /path/to/sampleB_R2.fastq.gz
sampleB PACBIO  /path/to/sampleB_hifi.fastq.gz
```

## Pre-assembled genomes

GABI also accepts pre-assembled genomes - in which case only limited QC data can be generated of course. 

A samplesheet for pre-assembled genomes looks as follows:

```TSV
sample  assembly
sampleA /path/to/sampleA.fasta
sampleB /path/to/sampleB.fasta
```

Then you can run GABI as usual:

=== "All default settings"

    ```BASH
    nextflow run bio-raum/gabi -profile my_profile \
    -r 1.4.0 \
    --input samples.tsv \
    --run_name MyAssemblies
    ```

=== "Skipping AMR predictions"

    Skip the optional AMR prediction steps

    ```BASH
    nextflow run bio-raum/gabi -profile my_profile \
    -r 1.4.0 \
    --input samples.tsv \
    --skip_amr \
    --run_name MyAssemblies
    ```