# Frequently asked questions

This section will be expanded as we find new questions and potential pitfalls. 

## General topics

### How is GABI different from other pipelines, like AQUAMIS?

GABI was heavily inspired by [AQUAMIS](https://gitlab.com/bfr_bioinformatics/AQUAMIS), so that question is quite sensible. GABI is fully implemented in Nextflow using a modular code design. With Nextflow come a range of "advantages" such as being able to translate the pipeline to various compute infrastructures and run with different kinds of software provisioning frameworks. It should also be significantly easier to extend than a pipeline written in Snakemake. Case in point, GABI supports not only Illumina short reads, but also Nanopore and Pacbio long reads - which means that you won't have to switch between different pipelines for different technologies. GABI also provides functionality beyond mere assembly, such as AMR profiling, serotyping and annotation. Finally, GABI is very strictly versioned and switching between versions does not require you to install separate versions of the pipeline - versioning is handled via Github and the nextflow `release` option (-r).

With all that said, GABI and AQUAMIS should behave fairly similarily on Illumina data as far as assembly and QC are concerned. 

### Is it OK to use Conda, or should I use a container manager?

**Please** do not use Conda if at all possible. We are making it an option in GABI because it has a broad user base. But Conda is *not* suitable for production purposes. 

Conda environments are prone to breaking; when you update your conda installation, or try to add additional packages. But more importantly, conda environments are not guaranteed to be reproducible. Many packages
have fairly lose dependencies. That's not condas' fault, but it means that a package you install today is not guaranteed to yield the same results as an installation performed a year from now - even with the exact same version! This is pretty much the opposite of "reproducibility" and a big no-no in scientific data analysis. 

Properly versioned containers on the other hand are fully reproducible - they are built once and can be re-used as many times and on as many systems as you like. They are always "the same". And your system does not have to perform time-consuming environment solving; you just download the container and are good to go.  

## Technologies

### Do any ONT reads work with GABI?

Well...maybe. But to be perfectly transparent, we are not testing with R9 reads or earlier and will not dedicate time to issues related to their support. With the latest R10 chemistry, results with GABI generally look good. We do recommend you use the super-accurate base calling (SUP) model 5.0.0 or later however, since our tests show this to be a prerequisite for consistently good assemblies. 

## Detecting contaminations

Generally speaking, detecting contaminations from read data is based on either the presence of variable sites when there shouldn't be any (bacteria are haploid, after all) or mixed taxonomic signals from e.g. Kmer analyses. Please see our more detailed discussion [here](contamination.md).