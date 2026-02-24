# Validation

This page describes the measures taken to check the validity of the GABI analysis, specifically for taxa with relevance to food safety (<i>Campylobacter</i>, <i>E. coli</i>, <i>L. monocytogenes</i> and <i>S. enterica</i>).

# Testing workflow reproducibility

It is common practice to verify the results of a workflow by repeating the analysis on the same data and using the same settings

- twice on the same system
- twice on a different system

A prefectly reproducible workflow will yield the exact same results across all runs. 

GABI has a built-in test that can be used for this purpose. 

```bash
nextflow run bio-raum/gabi -profile YOUR_PROFILE,test -r 1.0.1
```

where YOUR_PROFILE is set to whatever is appropriate for your system(s). This will automatically download a test data set (<i>E. coli</i>) and run the full workflow. You can verify the equivalence of the results in one of two ways:

## md5sum of the final assembly
```bash
md5sum results/samples/SAMEA2707761/assembly/SAMEA2707761.fasta
```
Expected result:
| Pipeline version | md5sum |
|------------------|--------|
| 1.0.1            | `e25cdf2cfc1e833ada0b29c9e5d6ad52` |

This approach will not cover the entirety of the analysis, but includes the "primary" result (i.e. the assembly). See the next option for when this can  nevertheless be preferrable. 


## md5sum of the final JSON
```bash
md5sum results/samples/SAMEA2707761/SAMEA2707761.qc.json
```

This check is comprehensive and recommended, but includes some caveats. Specifically, the final JSON includes, in addition to the various outputs, information about the executing user, the current date/time as well as the version of Nextflow used. If for some reason these are not identical between runs and systems, the md5sums will obviously not match. You can either remove that information from the JSON, or use the above mentioned assembly-level check.

## Additional caveats

Perfect reproducibility requires for all the relevant building blocks of the analysis to be "fix", i.e. identical. This is easily achieved using software containers (e.g. Apptainer, Singularity, Docker). It is however not guaranteed for Conda, where many packages have "loose" dependencies and can thus yield slightly different environments depending on Conda version and time of installation. If you find that when using Conda your md5sums do not match between your two systems, consider testing with a Container framework to make sure Conda isn't the culprit.

## Validation of the workflow

GABI has been validated for use with <i>Campylobacter</i>, <i>Escherichia coli</i>, <i>Listeria monocytogenes</i> and <i>Salmonella enterica</i> using publically available data.

### Against results from the AQUAMIS pipeline on ring trial data

GABI generates a range of results, but at its core it will perform a judgement call on whether a sample is safe to use, or whether it has any critical issues. There are, to our knowledge, no systematic benchmark data sets available to properly test this against a ground truth. Instead, GABI was run in tandem with the [AQUAMIS pipeline, v1.4.0 ](https://gitlab.com/bfr_bioinformatics/AQUAMIS) on extensive data from a published [ring trial](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2023.1253362/full) comprising over 500 Illumina libraries from 40 samples (10 per species). The status of each of these datsets was compared with the AQUAMIS pipeline to determine accuracy of the GABI results.

 Pipeline version | Sample status | Accuracy |
|------------------|--------| ----------- |
| 1.0.1            | [TSV](https://github.com/bio-raum/gabi/blob/main/assets/benchmark/gabi_vs_aquamis_1.0.1.tsv) | 0.97 |
| 1.2.0            | [TSV](https://github.com/bio-raum/gabi/blob/main/assets/benchmark/gabi_vs_aquamis_1.2.0.tsv) | 0.97 |

All 18 differences are related to the underlying algorithm for coverage calculation, where Aquamis uses an (invalid) approach that will count bases in overlapping reads [twice](https://bioinformatics.stackexchange.com/questions/5427/double-counting-coverage-of-overlapped-read-pairs) - which lifts the samples in question above the coverage thresholds and flags them as pass instead of fail.

### Against a contamination benchmark (Pightling et al, 2019)

GABIs ability to detect contaminations from sequence data has been validated for use with <i>Campylobacter</i>, <i>Escherichia coli</i>, <i>Listeria monoctogenes</i> and <i>Salmonella enterica</i> using data from a benchmark published in [Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1914-x), extended with [data](https://zenodo.org/records/4601406) published by the authors of AQUAMIS. 

#### Intraspecies contamination (cgMLST, ConfindR)


 Pipeline version | Species | Genetic distance |Recall | Precision |
|------------------|--------| ----------- | ---------- | --------- |
| 1.0.1, 1.2.0            | Campylobacter | >= 0.5 | 1.0 | 1.0
| 1.0.1, 1.2.0            | Campylobacter | 0.05 | 0.5 | 1.0
| 1.0.1, 1.2.0            | Campylobacter | 0.5 | 1.0 | 1.0
| 1.0.1, 1.2.0            | Campylobacter | 5 | 1.0 | 1.0
| 1.0.1, 1.2.0            | Escherichia coli | >= 0.5 | 1.0 | 1.0
| 1.0.1, 1.2.0            | Escherichia coli | 0.05 | 0.1 | 1.0
| 1.0.1, 1.2.0            | Escherichia coli | 0.5 | 1.0 | 1.0
| 1.0.1, 1.2.0            | Escherichia coli | 5 | 1.0 | 1.0
| 1.0.1, 1.2.0            | Listeria monocytogenes | >= 0.5 | 0.98 | 1.0
| 1.0.1, 1.2.0            | Listeria monocytogenes | 0.05 | 0.03 | 1.0
| 1.0.1, 1.2.0            | Listeria monocytogenes | 0.5 | 0.97 | 1.0
| 1.0.1, 1.2.0            | Listeria monocytogenes | 5 | 1.0 | 1.0
| 1.0.1, 1.2.0            | Salmonella enterica | >= 0.5 | 1.0 | 1.0
| 1.0.1, 1.2.0            | Salmonella enterica | 0.05 | 0.19 | 1.0
| 1.0.1, 1.2.0            | Salmonella enterica | 0.5 | 1.0 | 1.0
| 1.0.1, 1.2.0            | Salmonella enterica | 5 | 1.0 | 1.0

[Data - 1.0.1](https://github.com/bio-raum/gabi/blob/main/assets/benchmark/gabi_contamination_1.0.1.tsv), 
[Data - 1.2.0](https://github.com/bio-raum/gabi/blob/main/assets/benchmark/gabi_contamination_1.2.0.tsv)

These results suggest that intraspecies contaminants are robustly detectable at a genetic distance of 0.5% or higher. 

#### Interspecies contamination (cgMLST, ConfindR)

Pipeline version | Species | Contaminant | Recall | Precision |
| -------------- | ------- | ----------- | ------ | --------- |
| 1.0.1, 1.2.0          | Campylobacter | Escherichia coli | 1.0 | 1.0 |
| 1.0.1, 1.2.0          | Campylobacter | Listeria monocytogenes | 1.0 | 1.0 |
| 1.0.1, 1.2.0          | Campylobacter | Salmonella enterica | 1.0 | 1.0 |
| 1.0.1, 1.2.0          | Escherichia coli | Listeria monocytogenes | 1.0 | 1.0 |
| 1.0.1, 1.2.0          | Escherichia coli | Salmonella enterica | 1.0 | 1.0 |
| 1.0.1, 1.2.0          | Listeria monocytogenes | Escherichia coli | 1.0 | 1.0 |
| 1.0.1, 1.2.0          | Listeria monocytogenes | Salmonella enterica | 1.0 | 1.0 |
| 1.0.1, 1.2.0          | Salmonella enterica | Escherichia coli | 1.0 | 1.0 |
| 1.0.1, 1.2.0          | Salmonella enterica | Listeria monocytogenes | 1.0 | 1.0 |

\* Coverage of the contaminant is not displayed because all levels (>=10%) where correctly detected.