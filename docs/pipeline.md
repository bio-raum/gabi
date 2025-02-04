# Pipeline structure

![](../images/gabi_workflow.png)

Now, that looks very complex - and indeed it is. But we can break it down into a few key steps:

- Perform quality control of the read data, and merge libraries across lanes
- Group read data by sample id and check which assembly tool is appropriate based on the types of sequencing data we have available
- Assemble reads with the optimal tool
- Determine the species from the assembled genome sequence
- Perform quality checks on the assembly
- Perform genomic serotyping, if we have suitable tools available for that species
- Perform MLST typing on the assembly, if we have a pre-configured database for that species
- Annotate gene models in our assembled genome
- Predict antimicrobial resistance genes from our annotation
- Make a "pretty" QC report
