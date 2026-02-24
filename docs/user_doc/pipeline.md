# Pipeline structure

A DAG schema is included with this [code base][../dag-GABI.html]

Some of the key things that happen inside GABI:
- Perform quality control of the read data, and merge libraries across lanes
- Downsample reads to the global default (= 100X; can be adjusted)
- Group read data by sample id and check which assembly tool is appropriate based on the types of sequencing data we have available
- Assemble reads with the optimal tool (chain)
- Determine the species from the assembled genome sequence
- Perform quality checks on the assembly, including composition as contamination check
- Perform genomic serotyping, if we have suitable tools available for that species
- Perform MLST typing on the assembly, if we have a pre-configured database for that species
- Annotate gene models in our assembled genome
- Predict antimicrobial resistance genes from our assembly
- Call variants from the assembly (mostly for QC purposes, since there shouldn't be any!)
- Make a "pretty" QC report

# Does GABI distinguish between chromosomal assembly and "full" assembly, and what does that even mean?

Yes, GABI does make that distinction - meaning that some parts of the pipeline see the assembly with all the contigs, and some others only see the bits of the assembly that are likely to represent the bacterial chromosomes (without plasmids). 

| Pipeline section | Assembly used|
|------------------|--------------|
| AMR profiling    | All contigs  |
| MLST profiling   | All contigs  |
| Annotation       | All contigs  |
| Serotyping       | All contigs  |
| Reference genome | Chromosomes  |
| Variant calling  | Chromosomes  |
| Coverages        | All contigs  |
| Assembly QC      | All contigs  |

Why are we doing this? Basically, steps like idenfying the best reference genome match should not include plasmids, because these could throw off the algorithm. Likewise, if we want to know how complete our assembly is, this should ideally only include the chromosomes since the plasmids are generally "fluid" and not part of the "core" genome. That said, for most analysis steps, the entire assembly is used since the "biochemistry" of an isolate is determined by all the genes in the cell, not just those on the chromosomes. 
