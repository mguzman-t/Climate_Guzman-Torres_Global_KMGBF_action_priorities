#!/usr/bin/env python3
"""
Generate global KM-GBF policy layers with quantitative, score-selected primary
classes and subclasses.

The script generates all mapped interventions:

0   No eligible mapped intervention
1   Target 1: spatial planning
21  Target 2: additional connectivity-restoration need
22  Target 2: projected restoration
8   Target 8: climate-associated vulnerability
10  Target 10: production landscapes
31  Target 3: climate-resilient conservation in diffuse areas
32  Target 3: ecological OECM screening in intensified/channelized areas
33  Target 3: support existing PA/OECM coverage
34  Target 3: local area-based conservation in permeable areas
97  Unresolved exact score tie

How the primary class is selected
---------------------------------
1. Each class has an ecological eligibility mask.
2. Each eligible class receives a continuous raw score between 0 and 1.
3. Raw scores are converted to percentile scores from 0 to 100 separately
   within each class and scenario. This makes the class-specific scores
   comparable as RELATIVE PRIORITY, while preserving their different meanings.
4. In each grid cell, the eligible class with the highest percentile score is
   selected as the primary intervention.
5. If two or more eligible classes are tied within SCORE_TIE_TOLERANCE, the
   primary code is 97 rather than silently imposing a precedence hierarchy.

Important interpretation
------------------------
The score-selected primary class is the intervention for which the cell has the
highest relative priority compared with other cells eligible for that same
intervention. It is not a direct estimate of implementation feasibility, cost,
or probability of success.

OECM outputs are ecological screening products, not site-level eligibility or
designation results.
"""

import os
import glob
import json
import warnings

import numpy as np
import xarray as xr
import rioxarray


# =============================================================================
# CONFIGURATION
# =============================================================================

RICHNESS_DIR = (
    "/capacity/occr_davin/mguzman/P1/Connectivity/"
    "Lumb_sources_05grid"
)
CONNECTIVITY_INPUT_DIR = (
    "/capacity/occr_davin/mguzman/P1/Connectivity/"
    "results_omniscape_05grid/1_composite_agreement_maps"
)
HMOD_DIR = (
    "/capacity/occr_davin/mguzman/P1/Connectivity/hmod_05grid"
)
LU_ROOT = (
    "/capacity/occr_davin/mguzman/P1/ISIMIP3_LU_30yr_means/"
    "30yr_means"
)
PA_OECM_PATH = (
    "/capacity/occr_davin/mguzman/P1/Connectivity/"
    "WDPA_2026/pa_oecm_fraction_05deg.tif"
)
OUT_DIR = (
    "/capacity/occr_davin/mguzman/P1/Connectivity/"
    "kmgbf_05grid_quant_score"
)
os.makedirs(OUT_DIR, exist_ok=True)

SSPS = [
    "ssp126soc-adapt",
    "ssp370soc-adapt",
    "ssp585soc-adapt",
]
JUMPS = ["J1", "J2"]
LUMS = ["image", "magpie"]
RICHNESS_TAXA = ["amphibian", "bird", "reptile", "mammal"]
CONNECTIVITY_TAXA = ["bird", "mammal"]
BASELINE_PER = "1981_2010"
JUMP_TIME_MAP = {"J1": "2041_2070", "J2": "2071_2100"}

HIGH_VALUE_PERCENTILE = 66.0
BROAD_VALUE_PERCENTILE = 33.0
TOP_DECILE_PERCENTILE = 90.0
MIN_TOP_DECILE_TAXA = 2
MIN_VALID_RICHNESS_TAXA = 3
MIN_VALID_CHANGE_TAXA = 2
PA_OECM_GAP_THRESHOLD = 0.30
TAXON_DECLINE_THRESHOLD = -5.0
MAX_DECLINING_TAXON_FRACTION_TARGET3 = 0.25
MIN_DECLINING_TAXON_FRACTION_TARGET8 = 0.50
TARGET3_MEAN_CHANGE_STABILITY_THRESHOLD = -5.0
TARGET3_EXISTING_MEAN_DECLINE_THRESHOLD = -5.0

HMOD_LOW = 0.20
HMOD_HIGH = 0.55
LU_EPS = 0.001
SCORE_TIE_TOLERANCE = 1.0e-6
PRIMARY_TIE_CODE = 97

CRS = "EPSG:4326"
FLOAT_NODATA = -9999.0
UINT8_NODATA = 255

CLASS_CODES = [1, 21, 22, 8, 10, 31, 32, 34, 33]

PRIMARY_LABELS = {
    0: "No eligible mapped intervention",
    1: "T1 spatial planning",
    21: "T2 additional connectivity-restoration need",
    22: "T2 projected restoration",
    8: "T8 climate-associated vulnerability",
    10: "T10 production landscapes",
    31: "T3 climate-resilient conservation in diffuse areas",
    32: "T3 ecological OECM screening",
    33: "T3 support existing PA/OECM coverage",
    34: "T3 local area-based conservation",
    97: "Unresolved exact score tie",
}

MISSING_REASON_LABELS = {
    0: "Valid classification cell",
    1: "Missing PA/OECM coverage",
    2: "Missing connectivity",
    3: "Missing land-use transition",
    4: "Missing human modification",
    5: "Insufficient richness taxa",
    6: "Insufficient richness-change taxa",
}


# =============================================================================
# BASIC UTILITIES
# =============================================================================

def open_raster(path):
    return (
        rioxarray.open_rasterio(path, masked=True)
        .squeeze(drop=True)
        .drop_vars("band", errors="ignore")
    )


def align(data, template):
    return data.reindex_like(template, method="nearest")


def valid_mask(*arrays):
    if not arrays:
        raise ValueError("At least one DataArray is required.")
    output = arrays[0].notnull()
    for data in arrays[1:]:
        output = output & data.notnull()
    return output


def percentile_rank(data):
    """Return global percentile ranks from 0 to 100; ties use average ranks."""
    values = np.asarray(data.values, dtype=np.float64)
    valid = np.isfinite(values)
    output = np.full(values.shape, np.nan, dtype=np.float32)
    vector = values[valid]

    if vector.size == 0:
        return data.copy(data=output)
    if vector.size == 1:
        output[valid] = 100.0
        return data.copy(data=output)

    order = np.argsort(vector, kind="mergesort")
    sorted_values = vector[order]
    sorted_ranks = np.empty(vector.size, dtype=np.float64)

    start = 0
    while start < vector.size:
        end = start + 1
        while (
            end < vector.size
            and sorted_values[end] == sorted_values[start]
        ):
            end += 1
        sorted_ranks[start:end] = (start + end - 1) / 2.0
        start = end

    ranks = np.empty_like(sorted_ranks)
    ranks[order] = sorted_ranks
    output[valid] = (
        100.0 * ranks / (vector.size - 1)
    ).astype(np.float32)
    return data.copy(data=output)


def normalize_percentile(data, mask):
    return (
        percentile_rank(data.where(mask)) / 100.0
    ).clip(0.0, 1.0)


def write_float(data, path):
    output = (
        data.astype(np.float32)
        .rio.write_crs(CRS)
        .rio.write_nodata(FLOAT_NODATA)
    )
    output.rio.to_raster(
        path,
        dtype="float32",
        compress="LZW",
        nodata=FLOAT_NODATA,
    )


def write_uint8(data, path):
    output = xr.where(
        data.notnull(), data, UINT8_NODATA
    ).astype(np.uint8)
    output = (
        output.rio.write_crs(CRS)
        .rio.write_nodata(UINT8_NODATA)
    )
    output.rio.to_raster(
        path,
        dtype="uint8",
        compress="LZW",
        nodata=UINT8_NODATA,
    )


# =============================================================================
# RICHNESS
# =============================================================================

def match_richness_file(taxon, ssp, lum, period):
    folder = os.path.join(RICHNESS_DIR, taxon)
    if not os.path.isdir(folder):
        return None

    files = sorted(glob.glob(os.path.join(folder, "*.tif")))
    if period == BASELINE_PER:
        matches = [
            path for path in files
            if (
                "obsclim" in os.path.basename(path)
                or BASELINE_PER in os.path.basename(path)
            )
        ]
    else:
        matches = [
            path for path in files
            if (
                ssp in os.path.basename(path)
                and lum in os.path.basename(path)
                and period in os.path.basename(path)
            )
        ]

    if len(matches) > 1:
        warnings.warn(
            f"Multiple richness matches for {taxon}, {ssp}, {lum}, "
            f"{period}; using {matches[0]}"
        )
    return matches[0] if matches else None


def load_taxon_richness(taxon, ssp, lum, period, template=None):
    path = match_richness_file(taxon, ssp, lum, period)
    if path is None:
        print(f"  Missing richness: {taxon} | {ssp} | {lum} | {period}")
        return None
    data = open_raster(path).astype(np.float32)
    return align(data, template) if template is not None else data


def build_richness_evidence(ssp, lum, base_period, future_period):
    baseline = {}
    future = {}
    template = None

    for taxon in RICHNESS_TAXA:
        future_data = load_taxon_richness(
            taxon, ssp, lum, future_period, template
        )
        if future_data is None:
            return None
        if template is None:
            template = future_data
        future[taxon] = align(future_data, template)

        baseline_data = load_taxon_richness(
            taxon, ssp, lum, base_period, template
        )
        if baseline_data is None:
            return None
        baseline[taxon] = baseline_data

    percentiles = {
        taxon: percentile_rank(future[taxon])
        for taxon in RICHNESS_TAXA
    }
    percentile_stack = xr.concat(
        [percentiles[taxon] for taxon in RICHNESS_TAXA],
        dim="taxon",
    )
    valid_richness_taxon_count = (
        percentile_stack.notnull().sum("taxon").astype(np.float32)
    )
    sufficient_richness_data = (
        valid_richness_taxon_count >= MIN_VALID_RICHNESS_TAXA
    )
    vertebrate_value = (
        percentile_stack.mean("taxon", skipna=True)
        .where(sufficient_richness_data)
        .astype(np.float32)
    )
    top_decile_taxon_count = (
        (percentile_stack >= TOP_DECILE_PERCENTILE)
        .sum("taxon")
        .astype(np.float32)
    )

    changes = {}
    decline_flags = []
    change_valid_flags = []

    for taxon in RICHNESS_TAXA:
        change_valid = (
            baseline[taxon].notnull()
            & future[taxon].notnull()
            & (baseline[taxon] > 1.0)
        )
        change = xr.where(
            change_valid,
            100.0 * (future[taxon] - baseline[taxon]) / baseline[taxon],
            np.nan,
        ).astype(np.float32)
        changes[taxon] = change
        decline_flags.append(
            xr.where(
                change_valid,
                change < TAXON_DECLINE_THRESHOLD,
                0,
            ).astype(np.uint8)
        )
        change_valid_flags.append(change_valid.astype(np.uint8))

    decline_stack = xr.concat(decline_flags, dim="taxon")
    change_valid_stack = xr.concat(change_valid_flags, dim="taxon")
    change_stack = xr.concat(
        [changes[taxon] for taxon in RICHNESS_TAXA],
        dim="taxon",
    )

    valid_change_taxon_count = (
        change_valid_stack.sum("taxon").astype(np.float32)
    )
    sufficient_change_data = (
        valid_change_taxon_count >= MIN_VALID_CHANGE_TAXA
    )
    declining_taxon_count = (
        decline_stack.sum("taxon").astype(np.float32)
        .where(sufficient_change_data)
    )
    declining_taxon_fraction = xr.where(
        sufficient_change_data,
        declining_taxon_count / valid_change_taxon_count,
        np.nan,
    ).astype(np.float32)
    mean_change = (
        change_stack.mean("taxon", skipna=True)
        .where(sufficient_change_data)
        .astype(np.float32)
    )

    high_value = (
        (vertebrate_value >= HIGH_VALUE_PERCENTILE)
        | (top_decile_taxon_count >= MIN_TOP_DECILE_TAXA)
    ) & sufficient_richness_data
    broad_value = (
        vertebrate_value >= BROAD_VALUE_PERCENTILE
    ) & sufficient_richness_data
    stable_for_target3 = (
        sufficient_change_data
        & (
            declining_taxon_fraction
            <= MAX_DECLINING_TAXON_FRACTION_TARGET3
        )
        & (mean_change >= TARGET3_MEAN_CHANGE_STABILITY_THRESHOLD)
    )
    declining_for_target8 = (
        sufficient_change_data
        & (
            declining_taxon_fraction
            >= MIN_DECLINING_TAXON_FRACTION_TARGET8
        )
    )
    declining_richness = (
        sufficient_change_data
        & (declining_taxon_count >= 1)
    )
    mean_decline_existing = (
        sufficient_change_data
        & (mean_change < TARGET3_EXISTING_MEAN_DECLINE_THRESHOLD)
    )

    return {
        "template": template,
        "baseline": baseline,
        "future": future,
        "percentiles": percentiles,
        "changes": changes,
        "vertebrate_value_index": vertebrate_value,
        "top_decile_taxon_count": top_decile_taxon_count,
        "valid_richness_taxon_count": valid_richness_taxon_count,
        "sufficient_richness_data": sufficient_richness_data,
        "declining_taxon_count": declining_taxon_count,
        "declining_taxon_fraction": declining_taxon_fraction,
        "mean_taxon_richness_change_percent": mean_change,
        "valid_change_taxon_count": valid_change_taxon_count,
        "sufficient_change_data": sufficient_change_data,
        "high_biodiversity_value": high_value,
        "broad_biodiversity_value": broad_value,
        "stable_for_target3": stable_for_target3,
        "declining_for_target8": declining_for_target8,
        "declining_richness": declining_richness,
        "mean_decline_for_existing_coverage": mean_decline_existing,
    }


# =============================================================================
# COVERAGE AND CONNECTIVITY
# =============================================================================

def load_pa_oecm(template):
    if not os.path.exists(PA_OECM_PATH):
        raise FileNotFoundError(PA_OECM_PATH)
    return align(
        open_raster(PA_OECM_PATH).rio.write_crs(CRS), template
    ).astype(np.float32).clip(0.0, 1.0)


def load_hmod(ssp, jump, lum, template):
    path = os.path.join(
        HMOD_DIR,
        f"hmod_{jump}_{ssp}_{lum}_native05.tif",
    )
    if not os.path.exists(path):
        print(f"  Missing human modification: {path}")
        return None
    return align(open_raster(path).astype(np.float32), template)


def load_connectivity_slice(ssp, jump, lum):
    normalized_layers = []
    cumulative_layers = []
    template = None
    hmod_path = os.path.join(
        HMOD_DIR,
        f"hmod_{jump}_{ssp}_{lum}_native05.tif",
    )
    if not os.path.exists(hmod_path):
        return None
    hmod = open_raster(hmod_path).astype(np.float32)

    for taxon in CONNECTIVITY_TAXA:
        normalized_path = os.path.join(
            CONNECTIVITY_INPUT_DIR,
            f"consensus_mean_{taxon}_{ssp}_{lum}_{jump}_normcurrent.tif",
        )
        cumulative_path = os.path.join(
            CONNECTIVITY_INPUT_DIR,
            f"consensus_mean_{taxon}_{ssp}_{lum}_{jump}_cumcurrent.tif",
        )
        if (
            not os.path.exists(normalized_path)
            or not os.path.exists(cumulative_path)
        ):
            print(
                f"  Missing connectivity: {taxon} | {ssp} | {lum} | {jump}"
            )
            return None

        normalized = open_raster(normalized_path).astype(np.float32)
        cumulative = open_raster(cumulative_path).astype(np.float32)
        if template is None:
            template = normalized
            hmod = align(hmod, template)
        normalized_layers.append(align(normalized, template))
        cumulative_layers.append(align(cumulative, template))

    normalized_mean = xr.concat(
        normalized_layers, dim="connectivity_taxon"
    ).mean("connectivity_taxon", skipna=False)
    cumulative_mean = xr.concat(
        cumulative_layers, dim="connectivity_taxon"
    ).mean("connectivity_taxon", skipna=False)
    valid = valid_mask(normalized_mean, cumulative_mean, hmod)
    return (
        normalized_mean.where(valid),
        cumulative_mean.where(valid),
        hmod.where(valid),
    )


def derive_connectivity_thresholds():
    pooled_normalized = []
    pooled_cumulative = []
    for ssp in SSPS:
        for jump in JUMPS:
            for lum in LUMS:
                result = load_connectivity_slice(ssp, jump, lum)
                if result is None:
                    continue
                normalized, cumulative, hmod = result
                valid = valid_mask(normalized, cumulative, hmod).values
                pooled_normalized.append(normalized.values[valid])
                pooled_cumulative.append(cumulative.values[valid])

    if not pooled_normalized or not pooled_cumulative:
        raise FileNotFoundError("No complete connectivity inputs found.")

    normalized_values = np.concatenate(pooled_normalized)
    cumulative_values = np.concatenate(pooled_cumulative)
    thresholds = {
        "norm_q25": float(np.nanpercentile(normalized_values, 25)),
        "norm_q50": float(np.nanpercentile(normalized_values, 50)),
        "norm_q75": float(np.nanpercentile(normalized_values, 75)),
        "cum_q25": float(np.nanpercentile(cumulative_values, 25)),
    }
    with open(
        os.path.join(OUT_DIR, "connectivity_thresholds.json"),
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(thresholds, handle, indent=2, sort_keys=True)
    return thresholds


def classify_connectivity(ssp, jump, lum, thresholds, template):
    result = load_connectivity_slice(ssp, jump, lum)
    if result is None:
        return None
    normalized, cumulative, hmod = [
        align(data, template) for data in result
    ]
    valid = valid_mask(normalized, cumulative, hmod)
    low_flow = (cumulative < thresholds["cum_q25"]) & valid
    flow = (cumulative >= thresholds["cum_q25"]) & valid

    output = xr.zeros_like(normalized, dtype=np.uint8)
    output = xr.where(
        flow & (normalized > thresholds["norm_q75"]), 4, output
    )
    output = xr.where(
        flow
        & (normalized >= thresholds["norm_q50"])
        & (normalized <= thresholds["norm_q75"]),
        3,
        output,
    )
    output = xr.where(
        flow
        & (normalized >= thresholds["norm_q25"])
        & (normalized < thresholds["norm_q50"]),
        2,
        output,
    )
    output = xr.where(
        flow & (normalized < thresholds["norm_q25"]), 1, output
    )
    output = xr.where(
        low_flow & (hmod >= HMOD_HIGH), 7, output
    )
    output = xr.where(
        low_flow & (hmod >= HMOD_LOW) & (hmod < HMOD_HIGH),
        5,
        output,
    )
    output = xr.where(
        low_flow & (hmod < HMOD_LOW), 6, output
    )
    return output.where(valid)


# =============================================================================
# LAND USE
# =============================================================================

def load_lu(label, lum):
    if label == "histsoc":
        path = os.path.join(
            LU_ROOT,
            "histsoc",
            f"landuse_histsoc_mean_{BASELINE_PER}.nc",
        )
        if not os.path.exists(path):
            return None
        with xr.open_dataset(path) as dataset:
            return dataset.load()

    files = sorted(
        glob.glob(
            os.path.join(
                LU_ROOT,
                label,
                f"{label}_30y_means_{lum}_*.nc",
            )
        )
    )
    if not files:
        return None
    datasets = [xr.open_dataset(path).load() for path in files]
    try:
        return xr.concat(datasets, dim="member").mean(
            "member", skipna=True
        )
    finally:
        for dataset in datasets:
            dataset.close()


def land_classes(dataset, period):
    secondary = (
        dataset[f"secondary_forests_{period}"]
        + dataset[f"secondary_nonforests_{period}"]
    )
    return (
        secondary,
        dataset[f"pastures_{period}"],
        dataset[f"cropland_rainfed_{period}"],
        dataset[f"cropland_irrigated_{period}"],
        dataset[f"urbanareas_{period}"],
    )


def transition_evidence(base, future, base_period, future_period):
    (
        base_secondary,
        base_pasture,
        base_rainfed,
        base_irrigated,
        base_urban,
    ) = land_classes(base, base_period)
    (
        future_secondary,
        future_pasture,
        future_rainfed,
        future_irrigated,
        future_urban,
    ) = land_classes(future, future_period)

    rainfed = (future_rainfed - base_rainfed).clip(min=0)
    intensification = (
        future_irrigated - base_irrigated
    ).clip(min=0)
    pasture = (future_pasture - base_pasture).clip(min=0)
    urban = (future_urban - base_urban).clip(min=0)
    restoration = xr.where(
        (
            (future_rainfed - base_rainfed) < -LU_EPS
        )
        & (
            (future_secondary - base_secondary) > LU_EPS
        ),
        (future_secondary - base_secondary).clip(min=0),
        0.0,
    )

    total = rainfed + intensification + pasture + urban + restoration
    maximum = xr.concat(
        [urban, intensification, rainfed, pasture, restoration],
        dim="transition_type",
    ).max("transition_type", skipna=True)
    dominant = [
        (value == maximum) & (value > 0.5 * total) & (value > 0)
        for value in [
            urban,
            intensification,
            rainfed,
            pasture,
            restoration,
        ]
    ]
    stable = total <= LU_EPS

    category = xr.where(
        stable,
        5,
        xr.where(
            dominant[0],
            0,
            xr.where(
                dominant[1],
                1,
                xr.where(
                    dominant[2],
                    2,
                    xr.where(
                        dominant[3],
                        3,
                        xr.where(dominant[4], 4, 6),
                    ),
                ),
            ),
        ),
    ).squeeze(drop=True)

    conversion_magnitude = xr.concat(
        [urban, intensification, rainfed, pasture],
        dim="conversion_type",
    ).max("conversion_type", skipna=True)
    production_magnitude = xr.concat(
        [intensification, rainfed, pasture],
        dim="production_type",
    ).max("production_type", skipna=True)

    return {
        "category": category,
        "urban": urban.squeeze(drop=True),
        "intensification": intensification.squeeze(drop=True),
        "rainfed": rainfed.squeeze(drop=True),
        "pasture": pasture.squeeze(drop=True),
        "restoration": restoration.squeeze(drop=True),
        "conversion": conversion_magnitude.squeeze(drop=True),
        "production": production_magnitude.squeeze(drop=True),
    }


def get_landuse(ssp, lum, jump, template):
    historical = load_lu("histsoc", lum)
    future = load_lu(ssp, lum)
    if historical is None or future is None:
        return None

    if jump == "J1":
        evidence = transition_evidence(
            historical, future, BASELINE_PER, "2041_2070"
        )
    else:
        evidence = transition_evidence(
            future, future, "2041_2070", "2071_2100"
        )

    output = {}
    for name, data in evidence.items():
        rename = {
            old: new
            for old, new in [
                ("lat", "y"),
                ("lon", "x"),
                ("latitude", "y"),
                ("longitude", "x"),
            ]
            if old in data.dims
        }
        if rename:
            data = data.rename(rename)
        output[name] = align(data, template).astype(np.float32)
    return output


# =============================================================================
# ELIGIBILITY, SCORES, AND SCORE-SELECTED PRIMARY CLASS
# =============================================================================

def classify_and_score(
    richness,
    connectivity,
    landuse,
    coverage,
    valid,
):
    category = landuse["category"]

    high = richness["high_biodiversity_value"] & valid
    broad = richness["broad_biodiversity_value"] & valid
    stable = richness["stable_for_target3"] & valid
    declining8 = richness["declining_for_target8"] & valid
    declining = richness["declining_richness"] & valid
    mean_decline = (
        richness["mean_decline_for_existing_coverage"] & valid
    )

    gap = (coverage < PA_OECM_GAP_THRESHOLD) & valid
    covered = (coverage >= PA_OECM_GAP_THRESHOLD) & valid

    impeded = (connectivity == 1) & valid
    diffuse = (connectivity == 2) & valid
    intensified = (connectivity == 3) & valid
    channelized = (connectivity == 4) & valid
    limited = (connectivity == 5) & valid
    permeable = (connectivity == 6) & valid
    degraded = (connectivity == 7) & valid

    regional = intensified | channelized
    connectivity_restoration_need = impeded | limited | degraded
    stable_land_use = (category == 5) & valid
    restoration_transition = (category == 4) & valid
    production_change = category.isin([1, 2, 3]) & valid
    conversion_change = category.isin([0, 1, 2, 3]) & valid

    eligibility = {
        1: high & conversion_change,
        21: high & connectivity_restoration_need & ~restoration_transition,
        22: high & restoration_transition,
        8: high & declining8 & stable_land_use,
        10: broad & production_change,
        31: gap & high & stable & diffuse,
        32: gap & high & stable & regional,
        34: gap & broad & stable & permeable,
        33: covered & high & declining & mean_decline & conversion_change,
    }

    biodiversity = (
        richness["vertebrate_value_index"] / 100.0
    ).clip(0.0, 1.0)
    decline_breadth = richness[
        "declining_taxon_fraction"
    ].clip(0.0, 1.0)
    mean_change = richness[
        "mean_taxon_richness_change_percent"
    ]
    mean_decline_score = normalize_percentile(
        (-mean_change).clip(min=0.0), valid
    )
    mean_persistence_score = normalize_percentile(mean_change, valid)
    climate_vulnerability = (
        decline_breadth + mean_decline_score
    ) / 2.0
    climate_persistence = (
        (1.0 - decline_breadth).clip(0.0, 1.0)
        + mean_persistence_score
    ) / 2.0

    gap_score = (
        (PA_OECM_GAP_THRESHOLD - coverage)
        / PA_OECM_GAP_THRESHOLD
    ).clip(0.0, 1.0)

    connectivity_restoration_score = xr.where(
        degraded,
        1.00,
        xr.where(
            impeded,
            0.85,
            xr.where(limited, 0.70, 0.0),
        ),
    ).astype(np.float32)
    connectivity_conservation_score = xr.where(
        channelized,
        1.00,
        xr.where(
            intensified,
            0.85,
            xr.where(
                permeable,
                0.70,
                xr.where(diffuse, 0.60, 0.0),
            ),
        ),
    ).astype(np.float32)

    conversion_score = normalize_percentile(
        landuse["conversion"], conversion_change
    )
    production_score = normalize_percentile(
        landuse["production"], production_change
    )
    restoration_score = normalize_percentile(
        landuse["restoration"], restoration_transition
    )

    raw_scores = {
        1: (
            biodiversity
            + conversion_score
            + climate_vulnerability
        ) / 3.0,
        21: (
            biodiversity
            + connectivity_restoration_score
            + climate_vulnerability
        ) / 3.0,
        22: (
            biodiversity
            + restoration_score
            + connectivity_restoration_score
            + climate_vulnerability
        ) / 4.0,
        8: (
            biodiversity
            + decline_breadth
            + mean_decline_score
        ) / 3.0,
        10: (
            biodiversity
            + production_score
            + climate_vulnerability
        ) / 3.0,
        31: (
            biodiversity
            + climate_persistence
            + gap_score
        ) / 3.0,
        32: (
            biodiversity
            + connectivity_conservation_score
            + climate_persistence
            + gap_score
        ) / 4.0,
        34: (
            biodiversity
            + connectivity_conservation_score
            + climate_persistence
            + gap_score
        ) / 4.0,
        33: (
            biodiversity
            + climate_vulnerability
            + conversion_score
        ) / 3.0,
    }

    # Percentile-normalize each class among its eligible cells. These normalized
    # values are used for between-class selection in overlapping cells.
    normalized_scores = {
        code: normalize_percentile(
            raw_scores[code], eligibility[code]
        ).where(eligibility[code])
        for code in CLASS_CODES
    }

    score_stack = xr.concat(
        [normalized_scores[code] for code in CLASS_CODES],
        dim="intervention_code",
    ).assign_coords(intervention_code=CLASS_CODES)

    eligible_stack = xr.concat(
        [eligibility[code] for code in CLASS_CODES],
        dim="intervention_code",
    ).assign_coords(intervention_code=CLASS_CODES)

    eligible_count = eligible_stack.sum(
        "intervention_code"
    ).astype(np.float32)

    winning_score = score_stack.max(
        "intervention_code", skipna=True
    ).where(eligible_count > 0)

    winning_mask = (
        abs(score_stack - winning_score)
        <= SCORE_TIE_TOLERANCE
    ) & eligible_stack
    winning_count = winning_mask.sum(
        "intervention_code"
    ).astype(np.float32)

    # Use argmax only where the maximum is unique. Exact ties receive code 97.
    filled_stack = score_stack.fillna(-np.inf)
    winning_index = filled_stack.argmax("intervention_code")
    code_lookup = xr.DataArray(
        np.asarray(CLASS_CODES, dtype=np.uint8),
        dims=("intervention_code",),
        coords={"intervention_code": np.arange(len(CLASS_CODES))},
    )
    # Convert positional argmax indices to class codes.
    primary_values = np.take(
        np.asarray(CLASS_CODES, dtype=np.uint8),
        winning_index.values,
    )
    primary = winning_index.copy(
        data=primary_values
    ).astype(np.float32)
    primary = xr.where(eligible_count == 0, 0, primary)
    primary = xr.where(winning_count > 1, PRIMARY_TIE_CODE, primary)
    primary = primary.where(valid)

    sorted_scores = np.sort(
        filled_stack.values,
        axis=0,
    )
    best = sorted_scores[-1]
    second = sorted_scores[-2]
    best[np.isneginf(best)] = np.nan
    second[np.isneginf(second)] = np.nan
    runner_up_score = winning_score.copy(data=second.astype(np.float32))
    score_margin = winning_score - runner_up_score
    score_margin = xr.where(
        eligible_count == 1,
        winning_score,
        score_margin,
    ).astype(np.float32)

    target2_parent = eligibility[21] | eligibility[22]
    target3_parent = (
        eligibility[31]
        | eligibility[32]
        | eligibility[34]
        | eligibility[33]
    )
    target_overlap_count = (
        eligibility[1].astype(np.uint8)
        + target2_parent.astype(np.uint8)
        + target3_parent.astype(np.uint8)
        + eligibility[8].astype(np.uint8)
        + eligibility[10].astype(np.uint8)
    ).where(valid)

    output = {
        "primary_intervention": primary,
        "winning_normalized_score": winning_score,
        "runner_up_normalized_score": runner_up_score,
        "winning_score_margin": score_margin,
        "eligible_class_count": eligible_count.where(valid),
        "score_tie_count": winning_count.where(valid),
        "target_overlap_count": target_overlap_count,
        "biodiversity_score": biodiversity.where(valid),
        "climate_vulnerability_score": climate_vulnerability.where(valid),
        "climate_persistence_score": climate_persistence.where(valid),
        "coverage_gap_score": gap_score.where(valid),
        "connectivity_restoration_score": (
            connectivity_restoration_score.where(valid)
        ),
        "connectivity_conservation_score": (
            connectivity_conservation_score.where(valid)
        ),
        "conversion_magnitude_score": conversion_score.where(valid),
        "production_magnitude_score": production_score.where(valid),
        "projected_restoration_magnitude_score": (
            restoration_score.where(valid)
        ),
    }

    for code in CLASS_CODES:
        output[f"class_{code}_eligible"] = (
            eligibility[code].where(valid)
        )
        output[f"class_{code}_raw_score"] = (
            raw_scores[code].where(eligibility[code])
        )
        output[f"class_{code}_normalized_score"] = (
            normalized_scores[code]
        )

    return output


# =============================================================================
# SCENARIO BUILD AND OUTPUT
# =============================================================================

def build_missing_reason(coverage, connectivity, landuse, hmod, richness):
    reason = xr.zeros_like(coverage, dtype=np.uint8)
    reason = xr.where(coverage.isnull(), 1, reason)
    reason = xr.where(
        (reason == 0) & connectivity.isnull(), 2, reason
    )
    reason = xr.where(
        (reason == 0) & landuse.isnull(), 3, reason
    )
    reason = xr.where(
        (reason == 0) & hmod.isnull(), 4, reason
    )
    reason = xr.where(
        (reason == 0)
        & ~richness["sufficient_richness_data"],
        5,
        reason,
    )
    reason = xr.where(
        (reason == 0)
        & ~richness["sufficient_change_data"],
        6,
        reason,
    )
    return reason


def build_scenario(ssp, jump, lum, coverage_template, thresholds):
    base_period = (
        BASELINE_PER if jump == "J1" else "2041_2070"
    )
    future_period = JUMP_TIME_MAP[jump]
    richness = build_richness_evidence(
        ssp, lum, base_period, future_period
    )
    if richness is None:
        return None

    template = richness["template"]
    coverage = align(coverage_template, template)
    connectivity = classify_connectivity(
        ssp, jump, lum, thresholds, template
    )
    landuse = get_landuse(ssp, lum, jump, template)
    hmod = load_hmod(ssp, jump, lum, template)

    if connectivity is None or landuse is None or hmod is None:
        return None

    valid = valid_mask(
        coverage,
        connectivity,
        landuse["category"],
        hmod,
        richness["vertebrate_value_index"],
        richness["valid_change_taxon_count"],
    )
    valid = (
        valid
        & richness["sufficient_richness_data"]
        & richness["sufficient_change_data"]
    )

    classification = classify_and_score(
        richness,
        connectivity,
        landuse,
        coverage,
        valid,
    )
    missing_reason = build_missing_reason(
        coverage,
        connectivity,
        landuse["category"],
        hmod,
        richness,
    )

    dataset = xr.Dataset()
    for taxon in RICHNESS_TAXA:
        dataset[f"{taxon}_future_richness"] = richness["future"][taxon]
        dataset[f"{taxon}_richness_percentile"] = richness[
            "percentiles"
        ][taxon]
        dataset[f"{taxon}_richness_change_percent"] = richness[
            "changes"
        ][taxon]

    for name in [
        "vertebrate_value_index",
        "top_decile_taxon_count",
        "valid_richness_taxon_count",
        "declining_taxon_count",
        "declining_taxon_fraction",
        "mean_taxon_richness_change_percent",
        "valid_change_taxon_count",
        "sufficient_richness_data",
        "sufficient_change_data",
        "high_biodiversity_value",
        "broad_biodiversity_value",
        "stable_for_target3",
        "declining_for_target8",
        "declining_richness",
        "mean_decline_for_existing_coverage",
    ]:
        dataset[name] = richness[name].astype(np.float32)

    dataset["connectivity_class"] = connectivity.astype(np.float32)
    dataset["land_use_transition"] = landuse["category"]
    dataset["urban_expansion_magnitude"] = landuse["urban"]
    dataset["irrigated_intensification_magnitude"] = landuse[
        "intensification"
    ]
    dataset["rainfed_expansion_magnitude"] = landuse["rainfed"]
    dataset["pasture_expansion_magnitude"] = landuse["pasture"]
    dataset["projected_restoration_magnitude"] = landuse[
        "restoration"
    ]
    dataset["conversion_magnitude"] = landuse["conversion"]
    dataset["production_magnitude"] = landuse["production"]
    dataset["human_modification"] = hmod.astype(np.float32)
    dataset["pa_oecm_fraction"] = coverage.astype(np.float32)
    dataset["valid_data_coverage"] = xr.where(valid, 1.0, np.nan)
    dataset["missing_reason"] = missing_reason.astype(np.float32)

    for name, layer in classification.items():
        dataset[name] = layer.where(valid).astype(np.float32)

    dataset.attrs.update({
        "title": (
            "KM-GBF policy layers with score-selected primary intervention"
        ),
        "ssp": ssp,
        "land_use_model": lum,
        "time_jump": jump,
        "primary_labels": json.dumps(PRIMARY_LABELS),
        "primary_selection_method": (
            "Highest class-specific percentile score among ecologically "
            "eligible interventions"
        ),
        "tie_method": (
            "Code 97 when two or more eligible classes share the maximum "
            "normalized score within tolerance"
        ),
        "score_tie_tolerance": SCORE_TIE_TOLERANCE,
    })
    return dataset


def save_scenario(dataset, ssp, jump, lum):
    stem = f"{ssp}_{lum}_{jump}_native05"
    netcdf_path = os.path.join(
        OUT_DIR,
        f"kmgbf_score_selected_inputs_{ssp}_{lum}_{jump}.nc",
    )
    encoding = {
        variable: {
            "dtype": "float32",
            "_FillValue": FLOAT_NODATA,
            "zlib": True,
            "complevel": 4,
        }
        for variable in dataset.data_vars
    }
    dataset.to_netcdf(netcdf_path, encoding=encoding)

    write_uint8(
        dataset["primary_intervention"],
        os.path.join(
            OUT_DIR,
            f"kmgbf_primary_score_selected_{stem}.tif",
        ),
    )
    write_float(
        dataset["winning_normalized_score"],
        os.path.join(
            OUT_DIR,
            f"kmgbf_winning_normalized_score_{stem}.tif",
        ),
    )
    write_float(
        dataset["winning_score_margin"],
        os.path.join(
            OUT_DIR,
            f"kmgbf_winning_score_margin_{stem}.tif",
        ),
    )
    write_uint8(
        dataset["eligible_class_count"],
        os.path.join(
            OUT_DIR,
            f"kmgbf_eligible_class_count_{stem}.tif",
        ),
    )
    write_uint8(
        dataset["score_tie_count"],
        os.path.join(
            OUT_DIR,
            f"kmgbf_score_tie_count_{stem}.tif",
        ),
    )
    write_uint8(
        dataset["target_overlap_count"],
        os.path.join(
            OUT_DIR,
            f"kmgbf_target_overlap_count_{stem}.tif",
        ),
    )
    write_uint8(
        dataset["valid_data_coverage"],
        os.path.join(OUT_DIR, f"valid_data_coverage_{stem}.tif"),
    )
    write_uint8(
        dataset["missing_reason"],
        os.path.join(OUT_DIR, f"missing_reason_{stem}.tif"),
    )

    for code in CLASS_CODES:
        write_uint8(
            dataset[f"class_{code}_eligible"],
            os.path.join(
                OUT_DIR,
                f"class_{code}_eligible_{stem}.tif",
            ),
        )
        write_float(
            dataset[f"class_{code}_raw_score"],
            os.path.join(
                OUT_DIR,
                f"class_{code}_raw_score_{stem}.tif",
            ),
        )
        write_float(
            dataset[f"class_{code}_normalized_score"],
            os.path.join(
                OUT_DIR,
                f"class_{code}_normalized_score_{stem}.tif",
            ),
        )

    print(f"    Wrote {os.path.basename(netcdf_path)}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    print(
        "KM-GBF pipeline with score-selected primary interventions"
    )
    sample = load_taxon_richness(
        "bird", SSPS[0], LUMS[0], BASELINE_PER
    )
    if sample is None:
        raise FileNotFoundError(
            "Cannot build the spatial reference grid."
        )

    coverage = load_pa_oecm(sample)
    thresholds = derive_connectivity_thresholds()
    completed = 0

    for lum in LUMS:
        for ssp in SSPS:
            for jump in JUMPS:
                print(f"  {ssp} | {lum} | {jump}")
                dataset = build_scenario(
                    ssp,
                    jump,
                    lum,
                    coverage,
                    thresholds,
                )
                if dataset is None:
                    print("    skipped: missing required input")
                    continue
                save_scenario(dataset, ssp, jump, lum)
                dataset.close()
                completed += 1

    metadata = {
        "primary_labels": PRIMARY_LABELS,
        "class_codes": CLASS_CODES,
        "primary_selection_method": (
            "Highest class-specific percentile score among eligible classes"
        ),
        "normalization": (
            "Each class raw score is percentile-normalized among eligible "
            "cells in the same SSP-RCP, land-use model, and period"
        ),
        "tie_code": PRIMARY_TIE_CODE,
        "score_tie_tolerance": SCORE_TIE_TOLERANCE,
        "minimum_valid_richness_taxa": MIN_VALID_RICHNESS_TAXA,
        "minimum_valid_change_taxa": MIN_VALID_CHANGE_TAXA,
        "completed_scenario_combinations": completed,
    }
    with open(
        os.path.join(
            OUT_DIR,
            "kmgbf_score_selected_metadata.json",
        ),
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)

    print(f"Completed scenario combinations: {completed}")
    print(f"Outputs written to: {OUT_DIR}")


if __name__ == "__main__":
    main()
