# Y Chromosome Non-B DNA Annotation

Annotation of non-B DNA forming motifs across 141 chrY assemblies from T2T pangenome data (CEPH, HGSVC, v2 cohorts).

**Contact:** Marie Kratka  437364@gmail.com
**Location:** `/storage/praha5-elixir/projects/bioinf-fi/kratka/y_chr/`

---

## Where to Find the Annotations

All annotation outputs are in `data/`:

| Tool | Output location | Format | What it detects |
|------|----------------|--------|-----------------|
| g4Discovery | `data/g4Discovery/<sample>.bed` | BED | G-quadruplexes |
| non-B gfa (nBMST) | `data/nbmst/<sample>/` | GFF + TSV (14 files per sample) | G4, Z-DNA, IR, DR, MR, STR, A-phased repeats, cruciforms, triplexes, slipped DNA |
| G4Hunter | `data/g4hunter/Results_<sample>/<sample>.gff` | GFF | G-quadruplexes with propensity score |

Sample names correspond to accession IDs listed in `data/accessions.txt`.

---

## Input Data

141 chrY FASTA assemblies downloaded from the T2T pangenome S3 bucket:
- `T2T/scratch/chrY/ceph/v1/chrY_assemblies/`
- `T2T/scratch/chrY/hgsvc/v1/chrY_assemblies/`
- `T2T/scratch/chrY/v2/chrY_assemblies/`

Raw FASTAs: `data/accessions/<sample>.fa`

---

## Tools & Settings

### g4Discovery
- **Repo:** `workflow/tools/g4Discovery_on_metacentrum` (branch: `singularity-port`)
- **Script:** `workflow/scripts/2_run_g4discovery.sh`
- **Default settings**
- **Requires:** python environment, pqsfinder docker container from `/storage/praha5-elixir/projects/bioinf-fi/`

### non-B gfa (nBMST)
- **Repo:** `workflow/tools/non-B_gfa`
- **Script:** `workflow/scripts/3_run_nbmst.sh`
- **Settings:** Default settings, `-skipWGET`
- **Note:** Binary compiled on MetaCentrum, located at `workflow/tools/non-B_gfa/gfa`

### G4Hunter
- **Repo:** `workflow/tools/G4Hunter`
- **Script:** `workflow/scripts/4_run_g4hunter.sh`
- **Settings:** `-w 25 -s 1.2`
- **Requires:** `python27-modules` module on MetaCentrum

---

## Workflow

Steps 2-4 run as PBS job arrays on MetaCentrum (`-J 1-141`) with command `qsub <script.sh>`.

```
1_download_accessions.sh   # Download FASTAs from S3
2_run_g4discovery.sh       # G4 annotation via g4Discovery
3_run_nbmst.sh             # Broad non-B annotation via non-B gfa
4_run_g4hunter.sh          # G4 scoring via G4Hunter
```

Logs: `data/logs/<tool>_<array_index>.out/err`

---

## nBMST Output Files

Each sample in `data/nbmst/<sample>/` contains 14 files covering:

- `*_GQ.*` — G-quadruplexes
- `*_IR.*` — Inverted repeats (cruciform-forming)
- `*_DR.*` — Direct repeats (slipped DNA)
- `*_MR.*` — Mirror repeats (triplex-forming)
- `*_STR.*` — Short tandem repeats
- `*_Z.*` — Z-DNA
- `*_APR.*` — A-phased repeats (bent DNA)

Each motif type has both `.gff` and `.tsv` output.