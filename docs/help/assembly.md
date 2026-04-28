# Assembly strategies and options

## Disclaimer

GABI is an automated pipeline and while we try to incorporate state-of-the art processing strategies, they are not guaranteed to yield "perfect" results. But we suspect that GABI should get you *close enough* for many potential applications. Please perform your own analysis/validation if you plan on using GABI in a production scenario. 

## Long read assembly strategies

GABI has two alternative strategies for long read assembly - using either a single assembler step, or a multi-assembly pipeline and subsequent consensus finding. 

| Strategy | Assembler(s) | Command line option |
| -------- | ------------ | ------------------- |
| Single assembler | Flye |  |
| Multi-assembler | Canu, Flye, Miniasm, Necat, Raven, Plassembler | --autocycler |

To run a consensus assembly, GABI uses [Autocycler](https://github.com/rrwick/Autocycler) with a combination of [Canu](https://github.com/marbl/canu), [Flye](https://github.com/mikolmogorov/Flye), [Miniasm](https://github.com/lh3/miniasm), [Necat](https://github.com/xiaochuanle/NECAT), [Raven](https://github.com/lbcb-sci/raven) and [Plassembler](https://github.com/gbouras13/plassembler) - depending on the type of sequencing reads (ONT, Pacbio CLR or Pacbio HiFI) available.

Unsurprisingly, consensus assembly drastically increases run time and is recommended primarily if you want to get the absolutely best possible assembly out of GABI. You most likey also want to run this on a larger compute infrastructure (= high-performance cluster), unless you are really only trying to assembly one or two genomes. For many downstream applications, the single-tool approach may yield sufficiently accurate results - which is why it is the default option in GABI. We recommend you perform your own tests to see which strategy works best for your use case. 

## Assembly polishing

Polishing is a process by which an initial assembly draft is re-evaluated and improved by using a reference or raw read data. Depending on the sequencing technology, this may simply involve re-mapping the reads initially used for assembly to reconcile any dubious regions in the assembled sequence, or make use of complementary sequence data to fix errors that the assembled reads could not. A specific example would be the polishing of a long-read assembly with Illumina short reads to remove long-read specific issues such as homopolymer errors. 

| Data types | Polishing strategies |
| ---------- | -------------------- |
| Short reads | Polypolish (inside Shovill) |
| ONT | Medaka, Homopolish (optional) |
| ONT + short reads | Medaka, Polypolish |
| Pacbio | Homopolish (optional) |
| Pacbio + short reads | Polypolish |

## Technology-specific options

Generally, GABI runs fine with all-default settings. However, depending on your long read data, some adjustments may be necessary:

#### ONT data

`--porechop` Perform adapter trimming; this should not be necessary for more recently basecalled data. 

`--medaka_model`  The basecalling model used; only needed if your basecaller does not encode it in the sequence headers.

`--onthq` Use this option if your reads were basecalled with a SUP model

`--ont_min_q` Minimum quality ONT reads to keep

The last option is more meant to nudge the dataset towards "longer and better". Assemblers will perform downsampling of the reads anyway; but this  will normally not necessarily select for the "best" reads during that process. 

`--homopolish` Perform polishing of homopolymer errors using [Homopolish](https://github.com/ythuang0522/homopolish). Since this uses sequence information from related assemblies, some users may not wish to include such corrections.

#### Pacbio data

 `--pacbio_hifi` - use this option of your Pacbio data is from HiFi reads.

 `--homopolish` - Perform polishing of homopolymer errors using [Homopolish](https://github.com/ythuang0522/homopolish). Since this uses sequence information from related assemblies, some users may not wish to use such corrections. For Pacbio data, homopolish will only run when using CLR reads without short reads - i.e. is mutually exclusive with `--pacbio_hifi`.