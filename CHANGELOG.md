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
