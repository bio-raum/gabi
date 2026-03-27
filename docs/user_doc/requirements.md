# GABI Pipeline requirements

Below are some general guidelines to ensure that your data can be successfully analyzed by GABI. 

## Computing

Computing requirements are largely dependent on your input data. In principle, GABI should run on commodity hardware such as modern laptops or desktop computers. The minimum required specs include 4 CPU cores and 16 Gb of Ram (preferably > 32GB!) as well as ~50 GB of storage for the local reference data. That said, GABI was designed and optimized for powerful workstations and compute clusters, and such hardware is recommended if you plan on assembling many genomes - especially when using time-consuming consensus assembly for long-read data (`--autocycler`). 

While Nextflow, and consequently GABI, are technically compatible with Windows (through [WSL](https://learn.microsoft.com/en-us/windows/wsl/about)) and OSX, it is going to be easiest to run it on a Linux system. For more details, please see our [installation](installation.md) instructions. 

## Sequencing depth

GABI applies species-specific QC criteria to determine the suitabilty of a data set for analysis - one of which is the sequencing depth. As a rule of thumb, GABI expects somewhere between 20X to 40X of mean coverage (see [threshold](https://github.com/bio-raum/gabi/blob/main/assets/AQUAMIS_thresholds.json)), so it is recommended to aim for this range (ideally 50-100X) to ensure that your data does not trigger a warning or fail. 

## Contamination

GABI is very sensitive towards read contamination and will fail samples if it detects even lower levels of contamination. Please make sure you work with pure isolate cultures and that no contaminations are introduced during DNA extraction or library prep. Also note that some of the criteria applied to detect contamination may not perform optimally for Nanopore data, which is inherently more noisy. We are working to improve this. 

## Nanopore Reads

GABI supports processing of Nanopore (ONT) reads. Some recommendations include:

* Reads must be adapter-trimmed - and, if applicable, demultiplexed. GABI does not perform these processing steps.
* Basecalling should be performed with a recent version of [Dorado](https://github.com/nanoporetech/dorado) and a SUP (super-accurate) model
  * By default, GABI will use Medaka for polishing. This requires for your data to have been basecalled with Dorado. If this is not the case, use `--skip_medaka`
* If you have not yet concatenated the various individual FastQ files per sample, GABI can perform this task for you - just list one FastQ file per line in the sample sheet, each with the same sample ID.  

## IonTorrent Reads

GABI was neither designed nor tested to work with IonTorrent single-end reads. We simply do not see enough interest in this technology within our community. In principle, it should be possible to process Torrent data as if it were single-end Illumina data. However, proceed at your own risk as some tools may not appropriately deal with the torrent-specific error profiles. If you would like to see IonTorrent support for GABI, please open an issue on [github](https://github.com/bio-raum/gabi/issues). 

## Pacbio Reads

GABI supports processing of Pacbio reads. Some recommendations include:

* If at all possible, reads should be provided as HiFi (specify `--pacbio_hifi`) for the best performance and results
* GABI performs only rudimentary polishing of Pacbio assemblies as the best strategy is difficult to nail down(RS vs Sequel, HiFi vs subreads). Results overall should be good when using HiFi reads. 
* GABI does not perform any kind of trimming and demultiplexing of the data (use [Lima](https://lima.how/))

