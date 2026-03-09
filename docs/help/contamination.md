# Contamination

When sequencing bacterial isolates, avoiding contaminations is paramount. The introduction of "foreign" genetic material may notably degrade the quality of the assembly as well as introduce confounding signals that can result in incorrect conclusions across the various downstream analysis steps. For a discussion, see below.  

GABI uses several quality checks to ensure that the data and resulting assemblies are safe for used in downstream analyses. The underlying strategies are largely inspired by QC checks implemented in [Aquamis](https://gitlab.com/bfr_bioinformatics/AQUAMIS). 

## Intra-species contamination

Intra-species contamination may occur when, for example, cultures are produced from more than one bacterial clone - such as when picking overlapping colonies off a plate. Consequences of this kind of contamination can be hard to spot, especially when the original clones are closely related. The larger the genetic distance, the more likely users will notice e.g. an atypical genome size or fragmented assembies resulting from more complex assembly graphs. Methods such as MLST typing may yield conflicting or incorrect calls from a mix of genetic signatures. 

In addition to reference metrics for genome size and contig count, GABI employs [ConfindR](https://github.com/OLC-Bioinformatics/ConFindr) for select species to check for unexpected genetic polymorphisms in a set of broadly conserved reference genes as a proxy for intra-species contamination. 

## Inter-species contamination

Inter-species contamination may occur when a pure culture gets contaminated with one or several unrelated bacteria. One potential scenario is a failure to isolate one clonal isolate at the start, or insufficient care during sample handling. 

The consequences of inter-species contamination are typically more readily detectable, resulting in strong deviations from the expected assembly size and contig count, as well as affecting many other metrics such as GC content, gene count, etc. 

In addition to reference metrics for key assembly metrics, GABI performs taxonomic analysis of both raw reads as well as the final assembly to check whether the data indicates presence of more than one species. 

### Raw reads

Raw reads are analysed with Kraken2/Bracken to check for the presence of multiple species at an abundance of 5% or higher as an indicator for contamination. This metric is useful but not necessarily conclusive as e.g. highly abundant plasmids may trigger warnings when no actual contamination is present. 

#### Nanopore data

Nanopore data poses a particular challenge for the detection of contamination from SNP data since the read data is comparatively noisy. While GABI does try to perform such contamination checks on Nanopore data, the results are to be interpreted with a big grain of salt. Essentially, low levels of intra-species contaminations are unlikely to show up in Nanopore data since the (potentially) small number of genetic differences are drowned by the noise. To this end, we run ConfindR with rather stringent settings to prevent the noise from triggering warnings (a contaminating SNP must be supported by at least 5 reads, which would correspond to 10% when sequencing to a recommended depth of 50X; and only reads >= Q20 are used). Unfortunately, this still isn't a guarantee for a totally robust inference, depending on read depth and quality. In fact, as hinted at earlier, this strategy  will obscure cases of true contamination when read coverage or levels of contamination are low and converesely trigger contamination warnings for no reason when coverage is really high.

### Assembly

A more robust estimate of inter-species contamination may be found in the finished assembly. GABI checks each contig in the assembly using Kraken2 and computes the total number of nucleotides belonging to each identified best-hit species and genus. A contamination warning is trigger if more than 5% or the total assembly length is assigned to more than one species.  

Additionally, GABI also runs [CheckM2](https://github.com/chklovski/CheckM2), a dedicated tool for the detection of contamination from assemblies using machine learning. 