# 1.4.0

- Updated reference database to v1.2
- Added Autocycler to optionally run consensus assembly for long reads
- Moved redundant genome size calculations into its own sub workflow
- Defaulting genome size estimate to 6MB if calculation fails
- Medaka polishing: Check if reads contain basecalling model meta data; else skip if model isn't provided via command line
- Replaced ConfindR db download with proper installation workflow
- Replaced Chopper with Fastplong (faster, compatible with MultiQC)
- Removed read downsampling with RASUSA as all assembly tools/chains already take care of that
- Updated CheckM2 database download to specify a user-agent as a fix to aborting file transfers
- Adding pipeline settings to sample-level JSON
