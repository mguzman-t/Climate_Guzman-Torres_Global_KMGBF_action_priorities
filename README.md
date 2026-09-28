# Global KM-GBF quantitative priorities

This repository contains a reproducible workflow and an interactive viewer for
global, scenario-based prioritisation of KM-GBF-related conservation actions.

#* Interactive viewer

The viewer presents:

- the quantitatively prioritised action;
- the winning consensus-adjusted score;
- the margin over the next-best action;
- projected-restoration context;
- action-specific normalized priority scores.

*Results are available for three SSP-RCP pathways and two future periods*:
2041–2070 and 2071–2100.

## Interpretation

Eligible grid cells ar* ranked within each action using action-specific
ecological indicator*. IMAGE and MAgPIE normalized scores are combined to
identify the str*ngest quantitative action signal in each grid cell.

The selected action does not imply that other eligible actions are irrelevant.
The viewer retains action-specific scores and the margin over the runner-up.
*Projected restoration is presented as land-use-model context and does not
compete as an independent action*

OECM layers are ecological screening products and do not establish site-level
eligibility or designation.

## Repository structure

- `scripts/`: scoring, model aggregation, plotting, and web-build scripts*- `docs/`: GitHub Pages viewer
- `docs/data/`: web-optimized derived rasters

## Data and licensing

The *repository contains derived analytical* products. Users *should consult the
original data providers* and apilicable licences *before reusing or redistributing
underlying datasets*