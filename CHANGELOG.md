# 1.5.1

- Updating Medaka to 2.2.2 to add support for new Dorado basecalling model hac@6.0.0

# 1.5.0

- Increased phred threshold in ConfindR analysis of Illumina data from Q20 to Q25 to remove false positive contaminations reported to otherwise occur on current generaton Illumina instruments
- Removed Nanopore ConfindR as a fail category from final sample classification due to occurences of false contamination calls
- Combined serotyping results from ECTper and Abricate EcOH in the final HTML report as neither tool alone seems to do well on benchmark data
- Removed option to filter reads against host genome due to aggressive license restrictions of biobloomtools as well as unpredictable effects on downstream applications

# 1.4.2

- Added an upper version lock to the pipeline, preventing execution with new and untestet versions of nextflow
- Fixed an issue with Plassembler database installation
- Updated autocycler to 0.6.2

# 1.4.0

- Updated reference database to v1.2
- Added Autocycler to optionally run consensus assembly for long reads
- Added Plassembler to default long-read assembly workflows to improve recovery of plasmids when not using Autocycler
- Moved redundant genome size calculations into its own sub workflow
- Defaulting genome size estimate to 6MB if calculation fails
- Medaka polishing: Check if reads contain basecalling model meta data; else skip if model isn't provided via command line
- Replaced ConfindR db download with proper installation workflow
- Updated CheckM2 database download to specify a user-agent as a fix to aborting file transfers
- Adding pipeline settings to sample-level JSON
- Updated Chopper to version 0.12.0
