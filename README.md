# svt-cmv

Reproduction materials for the manuscript "Common method variance inverts the premise of stabilizer detection: Homogenization, noise-floor dominance, and the conservatism of sign-based aggregation" (under review). This repository holds the analysis code, Monte Carlo simulation outputs, and derived data tables required to regenerate every reported result.

## Contents

- `R/` - `analiz.R`, the single master script for the complete analysis pipeline.
- `python/` - `figures.py`, the figure-building script.
- `data/derived/` - derived result tables in CSV form: empirical cluster and moderator estimates, measurement-invariance stage profiles, latent panel results, threshold sensitivity checks, and the outputs of all Monte Carlo phases.
- `data/` - `svt_cmv_arsiv_public.rds`, the public archive object.
- `figures/` - Figures 1 to 3 in PNG, TIFF, and PDF formats.

## Reproduction notes

R 4.4 or newer is required, together with the lavaan, psych, boot, readxl, and stringi packages. The single master script `R/analiz.R` runs the empirical pipeline, all Monte Carlo phases, and the archive builder under master seed 9186. Tier-2 simulation cells are computationally heavy and take roughly 35 minutes each. Once the R pipeline has finished, `python/figures.py` builds Figures 1 to 3 at 300 dpi.

## Data statement

Participant-level survey data are not distributed in this repository. They are available from the corresponding author under the conditions stated in the manuscript. The derived tables provided here are sufficient to reproduce every table and figure reported in the paper.

## Licenses

Code in this repository is released under the MIT License. The derived data files are released under CC BY 4.0.

## Citation

Full citation information will be added when the paper is published.
