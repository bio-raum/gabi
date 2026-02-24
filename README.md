# bio-raum/gabi

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![run with apptainer](https://img.shields.io/badge/apptainer-run?logo=apptainer&logoColor=3EB049&label=run%20with&labelColor=000000)](https://apptainer.org/)

**G**enomic **A**nalysis of **B**acterial **I**solates, or GABI for short, is a bioinformatic workflow to assemble bacterial genomes, annotate genes, detect antimicrobial resistance genes, perform serotyping for select taxa as well assign MLS types. GABI supports both (Illumina) short reads and ONT/Pacbio long reads. Hybrid assembly approaches are supported for ONT/Illumina data sets and require for different read types to carry the same sample ID for grouping. Please see our [documentation](https://bio-raum.github.io/gabi/latest/user_doc/usage/) for detailed information.  

If you are really impatient, check out our [quickstart](https://bio-raum.github.io/gabi/latest/user_doc/quickstart/) guide. 

![schema](images/GABI_Schema_v11.png)

## Developer team

GABI is developed and supported by the following people and organizations:

[Marc Höppner](https://github.com/marchoeppner), Landeslabor Schleswig-Holstein, [LSH](https://www.schleswig-holstein.de/DE/landesregierung/ministerien-behoerden/LLABOR)

## Acknowledgements

We thank the developers of the [AQUAMIS](https://gitlab.com/bfr_bioinformatics/AQUAMIS) pipeline for making some of the building blocks on which GABI is based publically available - specifically the ConfindR database and validation data for <i>Campylobacter spp.</i> as well as the reference intervals for a broad range of bacteria to determine assembly status. In addition, we thank the [nf-co.re](https://nf-co.re/) community for developing the standards on which GABI is based. 