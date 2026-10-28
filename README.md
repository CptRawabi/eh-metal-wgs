# Genomic Characterization of a Multi-metal-Resistant *Enterobacter hormaechei*

## Overview

This project investigates the genomic characteristics of an environmental bacterial isolate, ZPb4, belonging to the species *Enterobacter hormaechei*.

The study combines whole-genome sequencing (WGS), phenotypic characterization, antimicrobial and metal resistance profiling, along with comparative genomic analyses to explore potential genetic determinants associated with metal tolerance and antibiotic resistance. The broader context is environmental AMR and the One Health interface between environmental microorganisms and human-associated bacterial populations. Genomic similarities are interpreted in an environmental and comparative context, without assuming direct transmission or establishing clinical risk from genomic evidence alone.

## Research objectives

- Characterize the draft genome and taxonomic placement of the environmental isolate.
- Investigate potential genetic determinants associated with lead, arsenic, and chromium tolerance.
- Identify and characterize genomic features associated with antimicrobial resistance.
- Compare the isolate with selected publicly available *Enterobacter* genomes.
- Assess plasmid replicons, insertion sequences, and other mobile genetic element-associated features.
- Document the bioinformatics workflows so that each step can be inspected and, where possible, repeated.

## Repository contents

```text
.
├── README.md
├── environments/              # one Conda environment file per analysis step
│   ├── eh_typing.yml
│   ├── prokka_env.yml
│   ├── blast_env.yml
│   ├── rgi_env.yml
│   ├── plasmidfinder.yml
│   ├── mge_env.yml
│   ├── ta_env.yml
│   ├── figures_env.yml
│   └── record/
│       └── base_original_state.yml   # historical record only, not for rebuilding
└── scripts/
    ├── original/              # command-line (WSL) scripts, kept exactly as written
    ├── figures/               # Python scripts for data organization and figures
    └── archive/               # scripts not used for the reported results
```

## Analysis steps, software and environments

Command-line steps were run in Ubuntu on WSL. Each `.yml` file in `environments/` recreates the software environment for one step.

| Step | Script (`scripts/original/`) | Environment | Software | Input / database |
|---|---|---|---|---|
| Genome collection and ST61 screen | `Genomic Data Collection.txt` | `eh_typing` | NCBI Datasets CLI 18.35.0 (conda; bulk download) and 18.36.0 (standalone binary, downloaded 22 Aug 2026; metadata and 28-genome download); jq 1.6; mlst 2.35.0 (`--scheme ecloacae`) | NCBI assemblies; PubMLST *E. cloacae* scheme bundled with mlst (snapshot 17 Jun 2026) |
| Genome annotation | `prokka.txt` | `prokka_env` | Prokka 1.15.6 (`--genus Enterobacter`), Prodigal 2.6.3, BLAST+ 2.16.0 (see limitations) | Assembled genomes (`.fna`) |
| Metal resistance search | `Bacmet.txt` | `blast_env` | BLAST+ 2.17.0 (`blastp`; identity ≥ 80% and query coverage ≥ 80%) | BacMet2 experimentally confirmed set (753 proteins); Prokka proteins |
| Antimicrobial resistance | `RGI_AMR.txt` | `rgi_env` | RGI 6.0.8, DIAMOND 2.2.5 (protein mode, default RGI cut-offs) | CARD 4.0.2 |
| Plasmid replicons | `PlasmidFinder.txt` | `plasmidfinder` | PlasmidFinder 2.1.6, KMA 1.6.13, BLAST+ 2.17.0 (`-d enterobacteriales`, minimum coverage 0.6, minimum identity 0.9) | `plasmidfinder_db` (git clone, commit dated 27 Nov 2025) |
| Mobile genetic elements | `Mobile Element Finder.txt` | `mge_env` | MobileElementFinder 1.1.2, Python 3.9.23 | Database bundled with the package (mgedb 1.1.1) |
| Toxin–antitoxin systems | `TASmania.txt` | `ta_env` | TAtouScan 0.2.0, HMMER 3.4 | Zenodo record 20059258, patched with `patch_db.sh` |
| Data organization and figures | `scripts/figures/` | `figures_env` | Python 3.13.13, pandas 3.0.5, numpy 2.5.2, matplotlib 3.11.1, seaborn 0.13.2, scipy 1.18.0, biopython 1.88, openpyxl 3.1.5 | Tables produced by the steps above |

### Steps run on Galaxy Europe

The following were run on the Galaxy Europe server (https://usegalaxy.eu) and are documented separately:

- Pairwise ANI with FastANI (v1.3, default parameters; species threshold ≥ 95%).
- Final MLST typing with MLST (v2.22.0, *E. cloacae* scheme, default parameters).
- Pan-genome analysis  using Panaroo (v1.8.0)
- The CDSs annotation with eggNOG-mapper14 (v2.1.13) against the eggNOGdatabase15 (v5.0.2), default parameters. 

Galaxy tool versions can change over time, so rerunning later may not give identical results.

## How to run the steps

1. Install Conda (Miniforge) in Ubuntu/WSL.
2. Create an environment from its file, for example:
   ```bash
   mamba env create -f environments/prokka_env.yml
   conda activate prokka_env
   ```
3. Set up the databases (see below).
4. Run the step's script with its matching environment active.

The scripts are archived exactly as they were run, so they contain machine-specific details. Work on a copy to re-run them:

- They use absolute paths such as `/mnt/d/...`, which must be adjusted.
- The text files carry Windows line endings. Paste commands into the terminal, or clean a copy first with `tr -d '\r' < script.txt > script.sh`.
- Activate the environment before running (the scripts contain no `conda activate` lines, except where stated in them).

### Databases

| Database | Source | Identifier |
|---|---|---|
| CARD | `rgi auto_load` inside `rgi_env` | v4.0.2. Checksums of the four FASTA index files: `card_database_v4.0.2.fasta` `500f594b42c6ac97fdbe17a2f788ffd8`; `card_database_v4.0.2_all.fasta` `aa8bb51771f095348d37c72b59b181e9`; `wildcard_database_v4.0.2.fasta` `e3a2d35e7cfb3e3a06a84e3972ca8a70`; `wildcard_database_v4.0.2_all.fasta` `3b126fe280ecc4cc5f8a845c6f970196` |
| BacMet2 (experimentally confirmed) | BacMet website (http://bacmet.biomedicine.gu.se) | 753 sequences, downloaded 10 Jul 2026, md5 begins `9d53ffbaf738428708a3121351101fe2 ` |
| PlasmidFinder database | `git clone` of `plasmidfinder_db` (Genomic Epidemiology, Bitbucket) | Commit `07648fddc4c6cf0330c9fc041dd6f52b4a69f8d6`, dated 27 Nov 2025 |
| TAtouScan database | Zenodo record 20059258 (`tatouscan_db.tar.gz`), then `patch_db.sh` | Folder fingerprint after patching: `9ec4cbfd62ef521ba1ddf0af5847ea6e`|
| MLST scheme | Bundled with mlst 2.35.0 (`ecloacae`) | Snapshot 17 Jun 2026, md5 begins `4ea718da80bc68561483b8317c37e316`|
| Reference genomes | NCBI, by accession | The 28 accessions are listed in `Genomic Data Collection.txt` |

## Manual and unscripted steps

Several decisions were made by hand and are not recorded in any script:

- **Strain selection.** All complete/chromosome-level *Enterobacter* assemblies were downloaded and typed by MLST to identify ST61. Strains relevant to metal tolerance and environmental origin were then selected manually in an Excel master sheet built from the metadata table. The 28 accessions in the script are the result of this selection. Of the 32 comparative genomes, 28 come from that list, and 4 (including 2 contig-level assemblies) were added manually. 
- **RGI summaries.** The Python summary block inside `RGI_AMR.txt` was not used for the reported results. The per-strain `rgi_*.txt` tables were concatenated and analysed in Excel.
- **PlasmidFinder summaries.** The combined Excel table and the coordinate table were created outside the saved scripts.
- **ZPb4.** The assembly of the study isolate is not downloaded by any script.

## Known limitations of the reproducibility record

1. Environment files were rebuilt to match the original conditions as far as the Conda history allowed. The original RGI, BLAST and Python work ran in a shared `base` environment that was later modified. `base_original_state.yml` records its current state for reference.
2. **Prokka.** The original run used Prokka 1.15.6. In the rebuilt environment, BLAST+ is 2.16.0 instead of 2.17.0 because the two versions could not be installed together with Prokka's Perl dependencies. Gene predictions come from Prodigal 2.6.3, which matches the original, and BLAST is expected to affect only gene product naming.
3. **BacMet.** BLAST+ 2.17.0 is inferred from the Conda history (the run took place before the `base` environment was changed), not from a recorded command.
4. **CARD.** Version 4.0.2 is taken from the original session's download. The database was reloaded on 25 Oct 2026 for verification. A verification re-run of RGI on regenerated proteins produced one extra hit in two strains and the discrepancy could not be verified as the database and tool versions were the same.
5. **TAtouScan database.** The released database lacks metadata rows for some profiles, which made the tool fail. `patch_db.sh` adds the missing rows with family `Unknown_Family`, and toxin or antitoxin type is assigned from the profile name prefix. 8 profiles were patched.
6. **MLST.** The command-line MLST run was used only to find ST61 candidates. The sequence types reported in the manuscript come from the Galaxy Europe run.
7. **Genome downloads.** NCBI queries return the database as it is on the day of the query, so a new run can return different genomes. The accession lists pin the dataset used.

## Data availability

Large sequencing files, intermediate files, and database resources are not stored in this repository. Public reference genomes are identified by NCBI accession number. The ZPb4 accession number will be given once the manuscript has been published.

## Project status

This repository documents an ongoing research project. Workflow validation, version documentation, and reproducibility checks are in progress, and not every step has been tested as part of an end-to-end pipeline.

## Acknowledgements

This work uses publicly available genomic resources and open-source bioinformatics software. The relevant tools, databases and publications are cited whereever applicable.
