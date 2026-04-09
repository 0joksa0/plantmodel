# Literature Harvest Checklist (Fit to Current Model)

## Goal
Collect literature datasets that can be directly compared to this model's outputs with minimal transformation, prioritizing:
- Arabidopsis thaliana
- photoperiod-resolved measurements
- sucrose/starch time-course (ZT)
- units compatible with `umol g^-1 FW` or convertible

## Fit Criteria (Scoring Guide)
- +30: Variable match (`sucrose`, `starch`, optionally `N`)
- +20: Time structure match (ZT/diel profile)
- +15: Photoperiod coverage (ideally 4/6/8/12h)
- +15: Unit compatibility (or clear conversion path)
- +10: Open, downloadable raw tables
- +10: Metadata quality (genotype, tissue, temperature, replicates)

## Current Inventory and Status
1. `data/paper/sugars_starch_reference.csv`
- Fit: `95/100`
- Status: ready
- Why: directly usable in current validation and plotting scripts.

2. `data/paper/N_data.csv`
- Fit: `65/100`
- Status: ready
- Why: useful for N supply-response validation, but not a diel time-course dataset.

3. `data/literature/raw/fairdom_fmv2_snapshot_inv123_s1.ro.zip`
- Fit: `85/100`
- Status: downloaded
- Why: broad FMv2 repository with experimental + simulation artifacts; strong candidate pool.

4. `data/literature/processed/fairdom_complete_data/lh2/...Metabolite_Data_(Golm_20-04-2015) (LH2)_YH.xlsx`
- Fit: `88/100`
- Status: extracted
- Why: likely direct metabolite curves, needs table normalization.

5. `data/literature/processed/fairdom_complete_data/lh3/...Metabolite_Data_(Golm_02-09-2015) (LH3).xlsx`
- Fit: `88/100`
- Status: extracted
- Why: likely direct metabolite curves, needs table normalization.

6. Dryad candidate: `doi:10.5061/dryad.6sf5184`
- Fit: `90/100`
- Status: blocked in automation (`401/403`)
- Why: appears highly relevant for sucrose under varied photoperiod/temperature conditions.

## Download Results (this session)
Downloaded:
- FAIRDOM snapshot archive (`investigation 123 snapshot 1`)
- Extracted LH2/LH3 complete-data ZIPs from that archive

Blocked:
- Dryad file streams (403/401 from non-interactive environment)
- FAIRDOM `investigation 74 snapshot 5` returned HTML (not archive) via direct non-interactive fetch

## Recommended Next Steps
1. Normalize FAIRDOM LH2/LH3 metabolite workbooks into CSV tables:
- target columns: `photoperiod_h`, `ZT`, `sucrose`, `starch`, `unit`, `tissue`, `genotype`

2. Add new normalized files under:
- `data/literature/processed/lh2_metabolites_normalized.csv`
- `data/literature/processed/lh3_metabolites_normalized.csv`

3. Extend plotting script to ingest multiple literature sources and show:
- per-source overlays
- pooled confidence envelope
- source-specific fit metrics

4. Manual fallback for Dryad:
- download files in browser session
- place them in `data/literature/raw/dryad_10_5061_dryad_6sf5184/`
- rerun normalization + overlay scripts.
