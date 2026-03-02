# GABI Pipeline requirements

Below are some general guidelines to ensure that your data can be successfully analyzed by GABI. 

## Computing

GABI is a comparatively light-weight pipeline and runs on a wide range of hardware, including your "average" laptop. At the same time, it is also capable of taking advantage of high-performance compute clusters, or "the cloud". 

At minimum, you will need 4 CPU cores, 16 GB of Ram and ~20.GB of disk space for the reference databases. In addition, the pipeline requires storage space for the intermediate files and the final results - which will depend on the size of your input data. 

While Nextflow, and consequently GABI, are technically compatible with Windows (through [WSL](https://learn.microsoft.com/en-us/windows/wsl/about)) and OSX, it is going to be easiest to run it on a Linux system. For more details, please see our [installation](installation.md) instructions. 

## Sequencing depth

GABI applies species-specific QC criteria to determine the suitabilty of a data set for analysis - one of which is the sequencing depth. As a rule of thumb, GABI expects somewhere between 20X to 40X of mean coverage (see [thresholds](https://github.com/bio-raum/gabi/blob/main/assets/AQUAMIS_thresholds.json)), so it is recommended to aim for this range to ensure that your data does not trigger a warning or fail. 

## Contamination

GABI is very sensitive towards read contamination and will fail samples if it detects even lower levels of contamination. Please make sure you work with pure isolate cultures and that no contaminations are introduced during DNA extraction or library prep. Also note that some of the criteria applied to detect contamination may not perform optimally for Nanopore data, which is inherently more noisy. We are working to improve this. 

## Nanopore Reads

GABI supports processing of Nanopore (ONT) reads. However, some restrictions apply.

* Reads should be generated with the R10 chemistry; we did not test nor recommend R9 reads for use with GABI.
* Reads should be of high quality (ideally Q20, generated with SUP basecalling); although the internally applied thresholds for filtering are user-configureable should your data be of lesser quality. 
* Reads must be adapter-trimmed (sequencing adapters, that is) - and, if applicable, demultiplexed. GABI does not perform these processing steps. 

## Pacbio Reads

At this point, GABI only supports Pacbio HiFi reads. This has been the standard output of PacBio sequencers for some time and helps to improve assembly quality and speed. If you still have legacy sub-read data, please see if you can convert them to CCS/HiFI reads with the appropriate [tool](https://ccs.how/).

## IonTorrent Reads

GABI was neither designed nor tested to work with IonTorrent single-end reads. We simply do not see enough interest in this technology within our community. In principle, it should be possible to process Torrent data as if it were single-end Illumina data. However, proceed at your own risk as some tools may not appropriately deal with the torrent-specific error profiles. If you would like to see IonTorrent support for GABI, please open an issue on [github](https://github.com/bio-raum/gabi/issues). 