# refs/ inventory (1376.72 GB, 46 lock-covered files)

| path | GB | files | category | rebuild |
|---|---|---|---|---|
| refs/age | 0.01 | 12 | age model artefacts | `scripts/train_*` |
| refs/age.bak_rclr | 0.01 | 5 | superseded age model artefacts (backup; not read) | `none` |
| refs/bacdive | 0.04 | 884 | BacDive API v2 response cache | `refilled on demand` |
| refs/cache | 0.24 | 541 | fast-load caches (regenerable) | `automatic` |
| refs/cmd | 0.02 | 64 | curatedMetagenomicData profiles + ExperimentHub index (scoring cohort) | `openbiota taxonomic-cohort` |
| refs/cohort | 9.88 | 736 | pathway reference cohort work dir (reads deleted after use) | `openbiota cohort` |
| refs/cohort_random | 146.46 | 664 | pathway reference cohort, uniform-random sampling work dir (retained reads; regenerable) | `openbiota cohort --sampling random` |
| refs/crosswalk | 0.03 | 9 | derived from GTDB/GlobDB/MetaPhlAn bridges | `make expansion-crosswalks` |
| refs/ctnpc | 0.0 | 5 | references | `make refs-ctnpc` |
| refs/db | 5.53 | 126 | built databases | `make build-db` |
| refs/decoys | 6.41 | 4 | food/background genomes | `make substrate-refs` |
| refs/expanded | 0.02 | 31 | locks + reconciliation (derived) | `make refs-expanded; make expansion` |
| refs/expansion_cohort | 0.0 | 103 | profiled locally from public reads; manifest inside each cohort file | `make expansion-cohorts` |
| refs/genomes | 6.71 | 17,915 | cache of reference genomes fetched on demand (NCBI/UHGG/GlobDB); each has a .source sidecar | `refilled on demand by openbiota.expansion.genomes` |
| refs/globdb_r232 | 584.09 | 346,262 | downloaded (locked); genomes/ unpacked from the archive | `make refs-expanded; make refs-globdb-genomes` |
| refs/gtdb | 0.49 | 13 | downloaded (locked) | `make refs-expanded` |
| refs/host | 6.52 | 8 | GRCh38 + PhiX bowtie2 host index | `openbiota build-host-index` |
| refs/kraken2 | 31.05 | 19 | downloaded (locked) publisher Kraken2/Bracken database | `make refs-expanded` |
| refs/metaphlan4_db | 91.4 | 217 | downloaded by MetaPhlAn (Jan26 + Jun23 indexes); marker dump and clade caches derived | `make refs-expanded (metaphlan --install); caches rebuild on first use` |
| refs/metaphlan_db | 3.45 | 11 | downloaded by MetaPhlAn 3 (scoring lane) | `make setup` |
| refs/micom | 12.19 | 3,166 | AGORA2 model library + community/LP caches | `make extension-refs` |
| refs/motus | 15.25 | 15 | downloaded (locked), unpacked | `make refs-expanded` |
| refs/mycobiome | 226.07 | 3,540 | mycobiome references | `make mycobiome-refs` |
| refs/nii | 0.0 | 0 | references | `make setup` |
| refs/panels | 0.02 | 23 | built DIAMOND panel databases | `make build-db` |
| refs/pathogens | 128.44 | 1,456 | pathogen target references | `make refs-strain` |
| refs/reference_cohort_samples.json | 0.0 | 1 | pathway cohort sample table | `openbiota cohort` |
| refs/reference_ranges.json | 0.0 | 1 | pathway reference ranges | `openbiota cohort` |
| refs/reference_ranges.prefix600k.json | 0.0 | 1 | superseded pathway ranges (prefix sampling; kept for audit) | `none` |
| refs/singlem | 21.56 | 793 | downloaded (locked) GlobDB metapackage | `make refs-expanded` |
| refs/strain_refs.lock.json | 0.0 | 1 | lock | `make strain-tools-lock` |
| refs/supplements | 19.62 | 7,008 | supplement references | `make substrate-refs` |
| refs/sylph | 51.33 | 19 | downloaded GTDB R232 sylph database | `make genome-lane` |
| refs/taxonomic_cohort.json | 0.01 | 1 | scoring cohort matrix (curatedMetagenomicData, MetaPhlAn 3) | `openbiota taxonomic-cohort` |
| refs/taxonomic_cohort_manifest.json | 0.0 | 1 | scoring cohort manifest | `openbiota taxonomic-cohort` |
| refs/validation | 9.88 | 111 | validation fixtures | `openbiota validate` |
