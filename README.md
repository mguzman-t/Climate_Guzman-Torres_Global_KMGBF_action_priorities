# Global KM-GBF strategic-action priorities

Interactive global viewer of spatial priorities for strategic actions contributing to the Kunming-Montreal Global Biodiversity Framework (KM-GBF).

The viewer summarises spatial prioritisation results across socioeconomic pathways (SSP-RCP), connectivity scales, species distribution models, climate models, land-use models, and connectivity-response assumptions.

## Interactive viewer

The web viewer allows comparison across:

- SSP1-2.6, SSP3-7.0 and SSP5-8.5;
- connectivity radii of 300, 700 and 1,000 km;
- prioritised strategic actions;
- ensemble agreement for individual actions;
- calibrated strategic-action scores;
- priority diagnostics;
- connectivity-class agreement; and
- restoration and recovery diagnostics.

## Ensemble design

Each SSP-RCP scenario and connectivity radius contains 180 ensemble members:

- 2 taxon-specific connectivity outputs;
- 3 species distribution models (SDMs);
- 5 global climate models (GCMs);
- 2 land-use models; and
- 3 resistance responses.

This gives:

```text
2 × 3 × 5 × 2 × 3 = 180 members
```

Agreement maps report the proportion of available ensemble members supporting a given result. Higher agreement therefore indicates greater consistency across the modelling ensemble.

## Strategic-action prioritisation

Eight strategic actions are considered:

- **T1**: spatial planning under biodiversity and intact-ecosystem loss;
- **T2**: ecosystem and connectivity restoration;
- **T8**: climate-adaptation for biodiversity conservation;
- **T10**: biodiversity-sensitive production landscapes;
- **T3**: protection of biodiversity-relevant and climate-persistent diffuse/permeable landscapes;
- **T3**: protection of landscapes relevant for ecosystem connectivity;
- **T3**: support for existing PA/OECM under climate vulnerability; and
- **T3**: safeguarding existing PA/OECM under land-use pressure.

The modal strategic-action map shows the action selected most frequently across available ensemble members in each grid cell.

Action-specific agreement maps show the percentage of available ensemble members for which an individual strategic action is selected as the highest-priority action.

## Score calibration

Strategic actions are compared using a fixed, positive-only, target-specific empirical calibration.

For each strategic action, positive scores are converted to empirical percentiles relative to all positive scores for the same action across SSPs, connectivity radii and ensemble members. To avoid bias, zeros are excluded from the empirical percentile reference distribution.

Calibration is performed once using the pooled reference data and the same calibration is subsequently applied across SSPs and connectivity radii.

The calibrated scores are used for member-level comparison among strategic actions. The original scores (geometric means) are retained in the analytical outputs but are not used directly for cross-action winner selection.

## Repository structure

```text
.
├── docs/
│   ├── index.html
│   ├── app.js
│   ├── catalog_png.json
│   ├── data_png/
│   └── data_boundaries/
├── scripts/
│   ├── config.py
│   ├── 01_validate_outputs.py
│   ├── 02_build_web_catalog.py
│   └── 03_render_png_layers.py
└── README.md
```

`docs/` contains the files required by the web viewer.

`docs/catalog.json` and the intermediate GeoTIFF copies in `docs/data/` are generated during the build process but are not needed to deploy the PNG-based viewer.

## Requirements

The processing scripts require Python 3 and the following Python packages:

```text
numpy
rasterio
Pillow
matplotlib
```

GDAL is optional. If `gdal_translate` is available, intermediate web GeoTIFFs are written as Cloud Optimised GeoTIFFs (COGs). Otherwise, the input GeoTIFFs are copied directly.

## Data availability

The processed PNG layers required for the interactive web viewer are included
in this repository.

The original analytical outputs underlying the KM-GBF prioritisation are not
included because of their size and are available from the authors upon
reasonable request.

## Validate the analytical outputs

Before rebuilding the viewer, validate the calibrated results:

```bash
python scripts/01_validate_outputs.py
```

The validation checks include the expected ensemble size, calibration metadata, raster grids and ranges, policy-selection percentages, modal priorities, agreement values, and consistency with member-level calibrated outputs.

## Build the viewer data

After validation, generate the web catalogue:

```bash
python scripts/02_build_web_catalog.py
```

Then render the PNG layers and final web catalogue:

```bash
python scripts/03_render_png_layers.py
```

The web application reads:

```text
docs/catalog_png.json
```

and the corresponding images under:

```text
docs/data_png/
```

## Run the viewer locally

From the repository root:

```bash
cd docs
python -m http.server 8000
```

Then open:

```text
http://localhost:8000/
```

## Spatial framework

Outputs are represented on the common 50-km Equal Earth analysis grid using the Equal Earth projection (EPSG:8857).

## Interpretation

The maps represent ensemble-based spatial priorities rather than deterministic predictions.

High ensemble agreement indicates that a result is supported consistently across the modelling ensemble. It does not imply complete certainty or establish that the selected strategic action is the only appropriate intervention for a grid cell.

Calibrated strategic-action scores represent the relative position of an action's positive scores within that action's pooled empirical distribution. They should not be interpreted as absolute percentages of ecological importance.

## Citation

Citation information will be added following publication of the associated study.