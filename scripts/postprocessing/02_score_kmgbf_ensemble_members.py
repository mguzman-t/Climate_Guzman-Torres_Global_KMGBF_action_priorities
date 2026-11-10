#!/usr/bin/env python3
"""
Create member-wise KM-GBF priority maps with direct path-driven inputs.

For each SSP x radius, the default run processes 180 matched members:

    2 connectivity taxa x 3 SDMs x 5 GCMs x 2 LUMs x 3 curves.

All inputs are resolved from the factor levels and path functions below.

Biodiversity horizon and redistribution
-----------------------------------------
Present and future richness are converted independently to taxon-balanced
spatial percentiles. Three biodiversity components are considered:

    B_current = present-day biodiversity importance
    B_future  = projected future biodiversity importance
    B_forward = max(B_current, B_future)

B_forward is an OR rule. It retains cells that are important
now or are projected to become important without allowing a large proportional
gain from a low baseline to dominate the biodiversity value variable.

Projected recipient areas are represented separately using richness change
calculated before any percentile transformation. For taxon t:

    DeltaR_t = Rfut_t - Rpres_t
    Gain_t   = max(DeltaR_t, 0)
    G_t      = clip(Gain_t / max(Rpres_t, 1), 0, 1)

Thus neither Gain_t nor G_t is a difference between future and current spatial
percentiles. Gain magnitude M_G is the mean G_t, and gain breadth Q_G is the
fraction of available taxa with Gain_t > 0. Recipient opportunity is:

    Recipient = GM(B_future, M_G, Q_G)

Two action diagnostics are also written:

    Recipient expansion   = GM(Recipient, PA/OECM coverage gap)
    Recipient restoration = GM(Recipient, HMOD) [impeded/diffuse]

These recipient diagnostics do not compete as an additional policy target.

Climate vulnerability
---------------------
For each taxon, proportional decline is:

    D = max((Rpres - Rfut) / Rpres, 0)

Mean decline magnitude M is kept on its absolute 0-1 scale (not spatially
reranked). Decline breadth Q is the fraction of available taxa with D > 0:

    V = sqrt(M * Q)
    persistence = 1 - V

Land-use common scaling
-----------------------
Positive land-use change magnitudes are divided by fixed pooled P99 references
read from common_scaling_references.json and clipped to [0, 1]. The reference
file must contain primary_loss, conversion_magnitude, production_magnitude and
secondary_recovery_proxy. Future production extent is already a 0-1 cell
fraction and is not P99-scaled.

Target scores
-------------
All components are on interpretable 0-1 scales and raw target scores are
geometric means:

    T1  = GM(B_current, primary_loss*, vulnerability)
    T2  = GM(B_forward, HMOD)                       [impeded/diffuse]
    T8  = GM(B_current, vulnerability)
    T10 = GM(B_forward, future_production_extent, production_expansion*)
    T31 = GM(B_future, intactness, coverage_gap, persistence) [diffuse]
    T32 = GM(B_future, coverage_gap)                [channelized/intensified]
    T33 = GM(B_current, coverage, vulnerability)
    T34 = GM(B_forward, coverage, conversion*)

Methodological notes
--------------------
Spatial percentiles express relative biodiversity importance within each
member, taxon and time horizon. Richness decline
and gain are therefore calculated from Rpres and Rfut before percentile
aggregation. Cells with Rpres = 0 are excluded from proportional magnitude
because a relative rate is undefined; their absolute gain still contributes to
gain breadth and remains available as a diagnostic. All competing targets use
the same availability mask and fixed
pooled P99 land-use references, so winner selection is on a common 0-100 scale.

Winner selection
----------------
Primary winners are selected directly from 100 * raw target scores. Final
scores are NOT spatially reranked. Separately reranked target hotspot scores
are produced only as within-target diagnostics and never affect winners.
Positive ties within TIE_TOLERANCE receive code 97; cells with no positive
available score receive code 0.

Outputs include member priority maps, member common-scale score archives,
modal priority and agreement, per-class selection percentages, mean common
scores, mean within-target hotspot percentiles, connectivity agreement,
biodiversity and recipient-opportunity diagnostics, Target-2 recovery
diagnostics, an inventory and metadata.
"""

import argparse
import csv
import itertools
import json
from pathlib import Path

import numpy as np
import rasterio

# =============================================================================
# PATHS
# =============================================================================

BASE = Path('/capacity/occr_davin/mguzman/P1/Connectivity/approach2')
RICHNESS_ROOT = BASE / 'gross_change_ee50km'
CONNECTIVITY_CLASS_ROOT = BASE / 'connectivity_agreement_sd4_ee50km/'
HMOD_ROOT = BASE / 'hmod_gcm_diagnostics'
LANDUSE_ROOT = BASE / 'landuse_transitions_ee50km'
SCALING_REFERENCE_PATH = LANDUSE_ROOT / 'common_scaling_references.json'
PA_OECM_PATH = (
    BASE / 'pa_oecm/pa_oecm_terrestrial_fraction_ee50km.tif'
)

OUTPUT_ROOT_ALL_TAXA = (
    BASE / 'kmgbf_memberwise_common_scale_180'
)

TEMPLATE_PATH = (
    RICHNESS_ROOT / 'bird/members/J1/ssp126'
    / 'bird_gam_gfdl-esm4_ssp126_J1_L_ee50km.tif'
)

# =============================================================================
# FACTOR LEVELS
# =============================================================================

SSPS = ['ssp126', 'ssp370', 'ssp585']
SSP_LANDUSE = {
    'ssp126': 'ssp126soc-adapt',
    'ssp370': 'ssp370soc-adapt',
    'ssp585': 'ssp585soc-adapt',
}
RADII = [6, 14, 20]
CONNECTIVITY_TAXA = ['bird', 'mammal']
SDMS = ['gam', 'glm', 'randomforest']
GCMS = [
    'gfdl-esm4',
    'ipsl-cm6a-lr',
    'mpi-esm1-2-hr',
    'mri-esm2-0',
    'ukesm1-0-ll',
]
LUMS = ['image', 'magpie']
CURVES = ['sensitive', 'linear', 'tolerant']
RICHNESS_TAXA = ['amphibian', 'bird', 'reptile', 'mammal']
JUMP = 'J1'

ARCHIVE_MEMBER_COUNT = (
    len(CONNECTIVITY_TAXA) * len(SDMS) * len(GCMS) * len(LUMS) * len(CURVES)
)

# =============================================================================
# OUTPUT CODES AND NUMERICAL SETTINGS
# =============================================================================

POLICY_CODES = [1, 21, 8, 10, 31, 32, 33, 34]
ALL_PRIMARY_CODES = [0, *POLICY_CODES, 97]
POLICY_LABELS = {
    0: 'No positive policy score',
    1: (
        'T1: biodiversity-inclusive spatial planning '
        'under land-use and climate co-exposure'
    ),
    21: 'T2: ecosystem and connectivity restoration priority',
    8: 'T8: climate-associated biodiversity vulnerability',
    10: 'T10: biodiversity-sensitive production landscapes',
    31: 'T3: maintain intact diffuse/permeable landscapes',
    32: 'T3: safeguard concentrated connectivity',
    33: (
        'T3: support existing PA/OECM coverage under '
        'climate-associated vulnerability'
    ),
    34: (
        'T3: safeguard existing PA/OECM coverage from '
        'projected land-use pressure'
    ),
    97: 'Unresolved exact score tie',
}
CONNECTIVITY_LABELS = {
    1: 'impeded',
    2: 'diffuse',
    3: 'channelized',
    4: 'intensified',
}

MIN_AVAILABLE_RICHNESS_TAXA = 2
MIN_AVAILABLE_CHANGE_TAXA = 2
TIE_TOLERANCE = 1.0e-6
NODATA_FLOAT = -9999.0
NODATA_UINT8 = 255

# =============================================================================
# GRID AND RASTER I/O
# =============================================================================


def load_template():
    if not TEMPLATE_PATH.is_file():
        raise FileNotFoundError(TEMPLATE_PATH)

    with rasterio.open(TEMPLATE_PATH) as source:
        values = (
            source.read(1, masked=True)
            .astype(np.float64)
            .filled(np.nan)
        )

        target = {
            'crs': source.crs,
            'transform': source.transform,
            'shape': source.shape,
            'width': source.width,
            'height': source.height,
            'profile': source.profile.copy(),
            'footprint': np.isfinite(values),
        }

    if target['crs'] != rasterio.crs.CRS.from_epsg(8857):
        raise ValueError('Template must use EPSG:8857')

    return target


def read_raster(path, target, minimum=None, maximum=None):
    if not path.is_file():
        raise FileNotFoundError(path)

    with rasterio.open(path) as source:
        if (
            source.count != 1
            or source.crs != target['crs']
            or source.transform != target['transform']
            or source.shape != target['shape']
        ):
            raise ValueError(f'Grid mismatch: {path}')

        values = (
            source.read(1, masked=True)
            .astype(np.float64)
            .filled(np.nan)
        )

    finite = np.isfinite(values)

    if not finite.any():
        raise ValueError(f'No finite values: {path}')

    if minimum is not None and np.any(values[finite] < minimum - 1e-7):
        raise ValueError(f'Values below {minimum}: {path}')

    if maximum is not None and np.any(values[finite] > maximum + 1e-7):
        raise ValueError(f'Values above {maximum}: {path}')

    return values


def write_float(path, array, target, description):
    profile = target['profile'].copy()
    profile.update(
        driver='GTiff',
        count=1,
        dtype='float32',
        nodata=NODATA_FLOAT,
        compress='LZW',
        tiled=True,
        blockxsize=256,
        blockysize=256,
    )

    temporary = path.with_name(path.stem + '.tmp.tif')

    with rasterio.open(temporary, 'w', **profile) as destination:
        destination.write(
            np.where(
                np.isfinite(array),
                array,
                NODATA_FLOAT,
            ).astype(np.float32),
            1,
        )
        destination.set_band_description(1, description)

    temporary.replace(path)


def write_uint8(path, array, target, description):
    profile = target['profile'].copy()
    profile.update(
        driver='GTiff',
        count=1,
        dtype='uint8',
        nodata=NODATA_UINT8,
        compress='LZW',
        tiled=True,
        blockxsize=256,
        blockysize=256,
    )

    temporary = path.with_name(path.stem + '.tmp.tif')

    with rasterio.open(temporary, 'w', **profile) as destination:
        destination.write(
            np.where(
                np.isfinite(array),
                array,
                NODATA_UINT8,
            ).astype(np.uint8),
            1,
        )
        destination.set_band_description(1, description)

    temporary.replace(path)

# =============================================================================
# CONTINUOUS TRANSFORMATIONS
# =============================================================================


def percentile_rank(array, mask=None):
    values = np.asarray(array, dtype=np.float64)
    valid = np.isfinite(values)

    if mask is not None:
        valid &= mask

    output = np.full(values.shape, np.nan, dtype=np.float64)
    vector = values[valid]

    if vector.size == 0:
        return output

    if vector.size == 1:
        output[valid] = 100.0
        return output

    order = np.argsort(vector, kind='mergesort')
    sorted_values = vector[order]
    sorted_ranks = np.empty(vector.size, dtype=np.float64)

    start = 0
    while start < vector.size:
        end = start + 1
        while end < vector.size and sorted_values[end] == sorted_values[start]:
            end += 1
        sorted_ranks[start:end] = (start + end - 1) / 2.0
        start = end

    ranks = np.empty_like(sorted_ranks)
    ranks[order] = sorted_ranks
    output[valid] = 100.0 * ranks / (vector.size - 1)
    return output


def percentile01(array, mask=None):
    return percentile_rank(array, mask=mask) / 100.0


def positive_percentile01(array, mask=None):
    """Rank positive values on 0-1 while preserving exact zeros as zero."""
    values = np.asarray(array, dtype=np.float64)
    valid = np.isfinite(values)
    if mask is not None:
        valid &= mask
    output = np.full(values.shape, np.nan, dtype=np.float64)
    output[valid & (values <= 0.0)] = 0.0
    positive = valid & (values > 0.0)
    if positive.any():
        ranked = percentile_rank(values, mask=positive) / 100.0
        output[positive] = ranked[positive]
    return output


def geometric_mean(arrays):
    stack = np.stack(arrays, axis=0)
    complete = np.all(np.isfinite(stack), axis=0)
    output = np.full(stack.shape[1:], np.nan, dtype=np.float64)

    clipped = np.clip(stack[:, complete], 0.0, 1.0)
    output[complete] = np.prod(clipped, axis=0) ** (1.0 / stack.shape[0])
    return output


def fixed_reference_scale(values, reference, mask=None):
    """Scale positive values by a fixed reference; preserve zero and nodata."""
    values = np.asarray(values, dtype=np.float64)
    available = np.isfinite(values)
    if mask is not None:
        available &= mask

    output = np.full(values.shape, np.nan, dtype=np.float64)
    output[available & (values <= 0.0)] = 0.0
    positive = available & (values > 0.0)
    output[positive] = np.clip(values[positive] / reference, 0.0, 1.0)
    return output


def load_scaling_references():
    """Read pooled P99 references required by the common-scale workflow."""
    if not SCALING_REFERENCE_PATH.is_file():
        raise FileNotFoundError(SCALING_REFERENCE_PATH)

    document = json.loads(SCALING_REFERENCE_PATH.read_text(encoding='utf-8'))
    variables = document.get('variables', document)
    required = (
        "primary_loss",
        'natural_loss',
        'conversion_magnitude',
        'production_magnitude',
        'secondary_recovery_proxy',
    )

    missing = [name for name in required if name not in variables]
    if missing:
        raise KeyError(
            f'Missing common-scaling references {missing} from '
            f'{SCALING_REFERENCE_PATH}. Rebuild the reference JSON with '
            'production_magnitude included before running this scorer.'
        )

    references = {}
    for name in required:
        record = variables[name]
        if not isinstance(record, dict) or 'p99' not in record:
            raise KeyError(
                f'p99 missing for {name} in {SCALING_REFERENCE_PATH}')
        reference = float(record['p99'])
        if not np.isfinite(reference) or reference <= 0.0:
            raise ValueError(f'Invalid P99 for {name}: {reference}')
        references[name] = reference

    return references, document

# =============================================================================
# MATCHED INPUT PATHS
# =============================================================================


def richness_path(taxon, sdm, gcm, ssp, metric):
    stem = f'{taxon}_{sdm}_{gcm}_{ssp}_{JUMP}_{metric}_ee50km.tif'
    return RICHNESS_ROOT / taxon / 'members' / JUMP / ssp / stem


def hmod_path(ssp, lum, gcm):
    return (
        HMOD_ROOT
        / f'hmod_landuse_J1_{SSP_LANDUSE[ssp]}_{lum}_{gcm}_ee50km.tif'
    )


def landuse_path(ssp, lum, gcm):
    return LANDUSE_ROOT / f'landuse_{ssp}_{lum}_{gcm}_J1.npz'


def connectivity_archive_path(ssp, radius):
    return (
        CONNECTIVITY_CLASS_ROOT
        / f'{ssp}_rad{radius}'
        / f'connectivity_classes_by_taxon_{ssp}_rad{radius}.npz'
    )

# =============================================================================
# RICHNESS: CONTINUOUS, TAXON-BALANCED
# =============================================================================


def build_richness_components(sdm, gcm, ssp, target):
    present_percentiles = []
    future_percentiles = []
    proportional_declines = []
    capped_proportional_gains = []
    positive_richness_gains = []

    for taxon in RICHNESS_TAXA:
        present = read_raster(
            richness_path(taxon, sdm, gcm, ssp, 'Rpres'),
            target,
            minimum=0.0,
        )
        future = read_raster(
            richness_path(taxon, sdm, gcm, ssp, 'Rfut'),
            target,
            minimum=0.0,
        )

        # Rank each taxon separately so taxa with larger richness ranges do not
        # dominate the multi-taxon biodiversity components.
        present_percentiles.append(
            percentile_rank(
                present,
                mask=target['footprint'] & np.isfinite(present),
            )
        )
        future_percentiles.append(
            percentile_rank(
                future,
                mask=target['footprint'] & np.isfinite(future),
            )
        )

        pair_available = (
            target['footprint']
            & np.isfinite(present)
            & np.isfinite(future)
        )

        # Calculate ecological change from the richness rasters themselves.
        # This is deliberately done before spatial percentile conversion.
        richness_change = np.full(
            target['shape'], np.nan, dtype=np.float64,
        )
        richness_change[pair_available] = (
            future[pair_available] - present[pair_available]
        )

        positive_richness_gain = np.full(
            target['shape'], np.nan, dtype=np.float64,
        )
        positive_richness_gain[pair_available] = np.maximum(
            richness_change[pair_available],
            0.0,
        )
        positive_richness_gains.append(positive_richness_gain)

        # Proportional decline requires a positive baseline. The absolute
        # 0-1 value is retained so SSP severity is not removed by reranking.
        proportional_change_available = pair_available & (present > 0.0)
        proportional_decline = np.full(
            target['shape'], np.nan, dtype=np.float64,
        )
        proportional_decline[proportional_change_available] = np.clip(
            -richness_change[proportional_change_available]
            / present[proportional_change_available],
            0.0,
            1.0,
        )
        proportional_declines.append(proportional_decline)

        # Gain magnitude uses positive DeltaR, divided by baseline richness to
        # place taxa with different pool sizes on a common 0-1 scale. The floor
        # and cap prevent very low baselines from producing unbounded leverage.
        capped_proportional_gain = np.full(
            target['shape'], np.nan, dtype=np.float64,
        )
        capped_proportional_gain[proportional_change_available] = np.clip(
            positive_richness_gain[proportional_change_available]
            / np.maximum(
                present[proportional_change_available],
                1.0,
            ),
            0.0,
            1.0,
        )
        capped_proportional_gains.append(capped_proportional_gain)

    def combine_taxa(percentile_arrays):
        stack = np.stack(percentile_arrays, axis=0)
        finite = np.isfinite(stack)
        count = np.sum(finite, axis=0)
        raw = np.divide(
            np.nansum(stack, axis=0),
            count,
            out=np.full(target['shape'], np.nan, dtype=np.float64),
            where=count >= MIN_AVAILABLE_RICHNESS_TAXA,
        )
        # Each taxon already contributes a 0-100 spatial percentile. Their
        # arithmetic mean is therefore the taxon-balanced biodiversity score;
        # a second spatial reranking would change that stated definition.
        score = raw / 100.0
        return score, raw, count.astype(np.float64)

    biodiversity_current, current_raw, current_count = combine_taxa(
        present_percentiles
    )
    biodiversity_future, future_raw, future_count = combine_taxa(
        future_percentiles
    )

    decline_stack = np.stack(proportional_declines, axis=0)
    decline_valid = np.isfinite(decline_stack)
    valid_change_count = np.sum(decline_valid, axis=0)

    vulnerability_raw = np.divide(
        np.nansum(decline_stack, axis=0),
        valid_change_count,
        out=np.full(target['shape'], np.nan, dtype=np.float64),
        where=valid_change_count >= MIN_AVAILABLE_CHANGE_TAXA,
    )

    # Absolute proportional decline preserves cross-SSP severity differences.
    decline_magnitude = np.clip(vulnerability_raw, 0.0, 1.0)
    decline_breadth = np.divide(
        np.sum((decline_stack > 0.0) & decline_valid, axis=0),
        valid_change_count,
        out=np.full(target['shape'], np.nan, dtype=np.float64),
        where=valid_change_count >= MIN_AVAILABLE_CHANGE_TAXA,
    )
    vulnerability = np.sqrt(decline_magnitude * decline_breadth)

    # Gain magnitude is the taxon mean of capped proportional gains. It is
    # derived from positive Rfut-Rpres, not from percentile differences.
    gain_stack = np.stack(capped_proportional_gains, axis=0)
    gain_available = np.isfinite(gain_stack)
    gain_magnitude_taxon_count = np.sum(gain_available, axis=0)
    gain_magnitude = np.divide(
        np.nansum(gain_stack, axis=0),
        gain_magnitude_taxon_count,
        out=np.full(target['shape'], np.nan, dtype=np.float64),
        where=gain_magnitude_taxon_count >= MIN_AVAILABLE_CHANGE_TAXA,
    )

    # Breadth is based directly on positive richness change and therefore also
    # recognises colonisation of cells with zero baseline richness. All paired
    # taxa form its denominator, including taxa with no gain.
    positive_gain_stack = np.stack(positive_richness_gains, axis=0)
    positive_gain_available = np.isfinite(positive_gain_stack)
    gain_breadth_taxon_count = np.sum(positive_gain_available, axis=0)
    gain_breadth = np.divide(
        np.sum(
            (positive_gain_stack > 0.0) & positive_gain_available,
            axis=0,
        ),
        gain_breadth_taxon_count,
        out=np.full(target['shape'], np.nan, dtype=np.float64),
        where=gain_breadth_taxon_count >= MIN_AVAILABLE_CHANGE_TAXA,
    )

    # This richness-unit diagnostic supports interpretation but does not enter
    # any common-scale policy score.
    gain_absolute_mean = np.divide(
        np.nansum(positive_gain_stack, axis=0),
        gain_breadth_taxon_count,
        out=np.full(target['shape'], np.nan, dtype=np.float64),
        where=gain_breadth_taxon_count >= MIN_AVAILABLE_CHANGE_TAXA,
    )

    # B_forward supports actions for which either present or future importance
    # is sufficient. Recipient opportunity remains a separate redistribution
    # diagnostic based on future importance and actual richness gains.
    biodiversity_forward = np.fmax(
        biodiversity_current,
        biodiversity_future,
    )
    recipient_opportunity = geometric_mean([
        biodiversity_future,
        gain_magnitude,
        gain_breadth,
    ])

    return {
        'biodiversity_current': biodiversity_current,
        'biodiversity_future': biodiversity_future,
        'biodiversity_forward': biodiversity_forward,
        'recipient_opportunity': recipient_opportunity,
        'gain_magnitude': gain_magnitude,
        'gain_breadth': gain_breadth,
        'gain_absolute_mean': gain_absolute_mean,
        'biodiversity_current_raw': current_raw,
        'biodiversity_future_raw': future_raw,
        'vulnerability': vulnerability,
        'vulnerability_raw': vulnerability_raw,
        'decline_magnitude': decline_magnitude,
        'decline_breadth': decline_breadth,
        'valid_current_richness_taxa': current_count,
        'valid_future_richness_taxa': future_count,
        'available_change_taxa': valid_change_count.astype(np.float64),
        'available_gain_magnitude_taxa': (
            gain_magnitude_taxon_count.astype(np.float64)
        ),
        'available_gain_breadth_taxa': (
            gain_breadth_taxon_count.astype(np.float64)
        ),
    }
# =============================================================================
# HMOD AND LAND-USE TRANSITIONS
# =============================================================================


def load_hmod(ssp, lum, gcm, target):
    return read_raster(
        hmod_path(ssp, lum, gcm),
        target,
        minimum=0.0,
        maximum=1.0,
    )


def load_landuse(ssp, lum, gcm, target):
    path = landuse_path(ssp, lum, gcm)
    if not path.is_file():
        raise FileNotFoundError(path)

    required = [
        "primary_loss",
        'natural_loss',
        'conversion_magnitude',
        'production_magnitude',
        'production_extent_future',
        'secondary_recovery_proxy',
    ]

    with np.load(path, allow_pickle=False) as archive:
        missing = [name for name in required if name not in archive.files]
        if missing:
            raise KeyError(f'Missing {missing} from {path}')
        output = {name: archive[name].astype(np.float64) for name in required}

    for name, values in output.items():
        if values.shape != target['shape']:
            raise ValueError(f'Shape mismatch for {name}: {path}')
        finite = np.isfinite(values)
        if np.any(values[finite] < -1e-7):
            raise ValueError(f'Negative {name}: {path}')
        if name == 'production_extent_future':
            maximum = float(
                values[finite].max()
            )

            # Fractions are expected to lie within [0, 1]. A small tolerance is
            # allowed for float32 round-off introduced when the aggregated
            # land-use arrays are written to the NPZ archives.
            if maximum > 1.0 + 1.0e-5:
                raise ValueError(
                    f'Materially invalid value for {name}: '
                    f'{path}; maximum = {maximum:.10f}'
                )

            output[name][finite] = np.clip(
                values[finite],
                0.0,
                1.0,
            )

    return output

# =============================================================================
# CONNECTIVITY CLASS ARCHIVE (TAXON-SPECIFIC, 180 MEMBERS)
# =============================================================================


def load_connectivity_classes(ssp, radius, target):
    path = connectivity_archive_path(ssp, radius)

    if not path.is_file():
        raise FileNotFoundError(
            f'Taxon-specific connectivity archive not found: {path}'
        )

    with np.load(path, allow_pickle=False) as archive:
        for name in ('connectivity_class', 'member_id', 'taxon'):
            if name not in archive.files:
                raise KeyError(f'{name} missing from {path}')

        classes = archive['connectivity_class'].astype(np.uint8)
        member_ids = archive['member_id'].astype(str)
        archive_taxa = archive['taxon'].astype(str)

    expected_shape = (ARCHIVE_MEMBER_COUNT, *target['shape'])

    if classes.shape != expected_shape:
        raise ValueError(
            f'Expected connectivity shape {expected_shape}, found {classes.shape}'
        )

    for taxon in CONNECTIVITY_TAXA:
        count = int(np.sum(archive_taxa == taxon))
        if count != ARCHIVE_MEMBER_COUNT // len(CONNECTIVITY_TAXA):
            raise ValueError(
                f'Unexpected member count for {taxon} in {path}: {count}'
            )

    valid_values = classes[classes != NODATA_UINT8]

    if not set(np.unique(valid_values)).issubset({1, 2, 3, 4}):
        raise ValueError(f'Unexpected connectivity class in {path}')

    lookup = {
        member_id: index
        for index, member_id in enumerate(member_ids)
    }

    if len(lookup) != len(member_ids):
        raise ValueError(f'Duplicate member identifiers in {path}')

    return classes, lookup, path


def member_identifier(taxon, sdm, gcm, lum, curve):
    return f'{taxon}__{sdm}__{gcm}__{lum}__{curve}'

# =============================================================================
# POLICY SCORING
# =============================================================================


def common_scale_policy_scores(
    connectivity_class,
    richness,
    hmod,
    landuse,
    coverage,
    references,
    target,
):
    # A shared mask makes the eight target scores directly comparable. A cell
    # is available only when every continuous input and the connectivity class
    # needed by the workflow are available.
    available = (
        target['footprint']
        & np.isfinite(richness['biodiversity_current'])
        & np.isfinite(richness['biodiversity_future'])
        & np.isfinite(richness['biodiversity_forward'])
        & np.isfinite(richness['vulnerability'])
        & np.isfinite(hmod)
        & np.isfinite(coverage)
        & np.isfinite(landuse['natural_loss'])
        & np.isfinite(landuse['conversion_magnitude'])
        & np.isfinite(landuse['production_magnitude'])
        & np.isfinite(landuse['production_extent_future'])
        & np.isfinite(landuse['secondary_recovery_proxy'])
        & np.isin(connectivity_class, [1, 2, 3, 4])
    )

    biodiversity_current = richness['biodiversity_current']
    biodiversity_future = richness['biodiversity_future']
    biodiversity_forward = richness['biodiversity_forward']
    recipient_opportunity = richness['recipient_opportunity']
    vulnerability = richness['vulnerability']
    persistence = np.clip(1.0 - vulnerability, 0.0, 1.0)
    intactness = np.clip(1.0 - hmod, 0.0, 1.0)
    coverage_gap = np.clip(1.0 - coverage, 0.0, 1.0)

    primary_loss = fixed_reference_scale(
        landuse["primary_loss"],
        references["primary_loss"],
        mask=available,
    )
    natural_loss = fixed_reference_scale(
        landuse['natural_loss'], references['natural_loss'], mask=available,
    )
    conversion = fixed_reference_scale(
        landuse['conversion_magnitude'],
        references['conversion_magnitude'],
        mask=available,
    )
    production_expansion = fixed_reference_scale(
        landuse['production_magnitude'],
        references['production_magnitude'],
        mask=available,
    )
    recovery = fixed_reference_scale(
        landuse['secondary_recovery_proxy'],
        references['secondary_recovery_proxy'],
        mask=available,
    )
    production_extent_future = np.full(
        target['shape'], np.nan, dtype=np.float64,
    )
    production_extent_future[available] = np.clip(
        landuse['production_extent_future'][available], 0.0, 1.0,
    )

    raw_scores = {
        # T1: current biodiversity exposed to simultaneous projected
        # natural-land contraction and richness decline.
        1: geometric_mean([
            biodiversity_current,
            primary_loss,
            vulnerability,
        ]),

        # T2: restoration need where biodiversity is important now or is
        # projected to become important.
        21: geometric_mean([
            biodiversity_forward,
            hmod,
        ]),

        # T8: current biodiversity exposed to projected climate-associated
        # richness decline.
        8: geometric_mean([
            biodiversity_current,
            vulnerability,
        ]),

        # T10: biodiversity-sensitive management where future production
        # landscapes and agricultural expansion coincide.
        10: geometric_mean([
            biodiversity_forward,
            production_extent_future,
            production_expansion,
        ]),

        # T31: intact, persistent and under-covered diffuse landscapes
        # expected to support future biodiversity.
        31: geometric_mean([
            biodiversity_future,
            intactness,
            coverage_gap,
            persistence,
        ]),

        # T32: concentrated connectivity in under-covered future-biodiversity
        # areas.
        32: geometric_mean([
            biodiversity_future,
            coverage_gap,
        ]),

        # T33: existing PA/OECM coverage containing current biodiversity
        # exposed to climate-associated decline.
        33: geometric_mean([
            biodiversity_current,
            coverage,
            vulnerability,
        ]),

        # T34: existing PA/OECM coverage where biodiversity important now or
        # in the future is exposed to projected gross conversion pressure.
        34: geometric_mean([
            biodiversity_forward,
            coverage,
            conversion,
        ]),
    }

    # Connectivity is the only categorical ecological gate.
    raw_scores[21][~np.isin(connectivity_class, [1, 2])] = np.nan
    raw_scores[31][connectivity_class != 2] = np.nan
    raw_scores[32][~np.isin(connectivity_class, [3, 4])] = np.nan
    for code in POLICY_CODES:
        raw_scores[code][~available] = np.nan

    # Common scores determine the winner. There is no final target reranking.
    common_scores = {code: 100.0 * raw_scores[code] for code in POLICY_CODES}

    # Within-target spatial percentiles are diagnostics only.
    hotspot_scores = {
        code: 100.0 * positive_percentile01(
            raw_scores[code], mask=np.isfinite(raw_scores[code]),
        )
        for code in POLICY_CODES
    }

    score_stack = np.stack([common_scores[code] for code in POLICY_CODES])
    filled = np.where(np.isfinite(score_stack), score_stack, -np.inf)
    best_index = np.argmax(filled, axis=0)
    best_score = np.max(filled, axis=0)
    any_score = np.isfinite(score_stack).any(axis=0)

    # A common-scale margin identifies weak or structurally dominant winners
    # without changing the selected class.
    second_best = np.partition(filled, -2, axis=0)[-2]
    winning_margin = np.full(target['shape'], np.nan, dtype=np.float64)
    margin_available = np.isfinite(best_score) & np.isfinite(second_best)
    winning_margin[margin_available] = (
        best_score[margin_available] - second_best[margin_available]
    )

    winner_count = np.sum(
        np.isfinite(score_stack)
        & np.isclose(
            score_stack,
            best_score[None, :, :],
            rtol=0.0,
            atol=TIE_TOLERANCE,
        ),
        axis=0,
    )

    winner = np.full(target['shape'], np.nan, dtype=np.float64)
    positive = any_score & (best_score > TIE_TOLERANCE)
    winner[available & ~positive] = 0
    unique = available & positive & (winner_count == 1)
    winner[unique] = np.take(
        np.asarray(POLICY_CODES, dtype=np.float64), best_index[unique],
    )
    winner[available & positive & (winner_count > 1)] = 97

    # Redistribution diagnostics are deliberately kept outside the policy
    # competition. They identify complementary opportunities for future
    # PA/OECM expansion and restoration in projected richness-recipient areas.
    recipient_expansion_opportunity = geometric_mean([
        recipient_opportunity,
        coverage_gap,
    ])
    recipient_restoration_opportunity = geometric_mean([
        recipient_opportunity,
        hmod,
    ])

    # Recipient diagnostics have their own availability condition and do not
    # suppress otherwise comparable policy scores where gain data are absent.
    recipient_available = available & np.isfinite(recipient_opportunity)
    recipient_masked = recipient_opportunity.copy()
    recipient_masked[~recipient_available] = np.nan
    recipient_expansion_opportunity[~recipient_available] = np.nan
    recipient_restoration_opportunity[
        (~recipient_available) | (~np.isin(connectivity_class, [1, 2]))
    ] = np.nan

    recovery_alignment = common_scores[21] * recovery
    restoration_gap = common_scores[21] * (1.0 - recovery)

    return {
        'available': available,
        'winner': winner,
        'winning_margin': winning_margin,
        'common_scores': common_scores,
        'hotspot_scores': hotspot_scores,
        'recipient_opportunity': 100.0 * recipient_masked,
        'recipient_expansion_opportunity': (
            100.0 * recipient_expansion_opportunity
        ),
        'recipient_restoration_opportunity': (
            100.0 * recipient_restoration_opportunity
        ),
        'recovery': recovery,
        'recovery_alignment': recovery_alignment,
        'restoration_gap': restoration_gap,
    }
# =============================================================================
# SSP x RADIUS PROCESSING
# =============================================================================


def process_group(ssp, radius, target, coverage, references, taxa, output_root):
    connectivity_stack, connectivity_lookup, connectivity_source = (
        load_connectivity_classes(ssp, radius, target)
    )

    members = list(itertools.product(taxa, SDMS, GCMS, LUMS, CURVES))
    n_members = len(members)

    richness_cache = {}
    hmod_cache = {}
    landuse_cache = {}

    winners = []
    available_members = []
    connectivity_members = []
    score_members = {code: [] for code in POLICY_CODES}
    hotspot_sums = {
        code: np.zeros(target['shape'], dtype=np.float64)
        for code in POLICY_CODES
    }
    hotspot_counts = {
        code: np.zeros(target['shape'], dtype=np.uint16)
        for code in POLICY_CODES
    }
    biodiversity_names = (
        'biodiversity_current',
        'biodiversity_future',
        'biodiversity_forward',
        'recipient_opportunity',
        'gain_magnitude',
        'gain_breadth',
    )
    biodiversity_sums = {
        name: np.zeros(target['shape'], dtype=np.float64)
        for name in biodiversity_names
    }
    biodiversity_counts = {
        name: np.zeros(target['shape'], dtype=np.uint16)
        for name in biodiversity_names
    }
    opportunity_names = (
        'recipient_expansion_opportunity',
        'recipient_restoration_opportunity',
    )
    opportunity_sums = {
        name: np.zeros(target['shape'], dtype=np.float64)
        for name in opportunity_names
    }
    opportunity_counts = {
        name: np.zeros(target['shape'], dtype=np.uint16)
        for name in opportunity_names
    }
    absolute_gain_sum = np.zeros(target['shape'], dtype=np.float64)
    absolute_gain_count = np.zeros(target['shape'], dtype=np.uint16)
    margin_sum = np.zeros(target['shape'], dtype=np.float64)
    margin_count = np.zeros(target['shape'], dtype=np.uint16)
    alignment_members = []
    gap_members = []
    recovery_members = []
    inventory = []

    group_dir = output_root / f'{ssp}_rad{radius}'
    member_output_dir = group_dir / 'member_priority_maps'
    member_output_dir.mkdir(parents=True, exist_ok=True)

    for member_index, (taxon, sdm, gcm, lum, curve) in enumerate(members):
        identifier = member_identifier(taxon, sdm, gcm, lum, curve)

        if identifier not in connectivity_lookup:
            raise KeyError(
                f'{identifier} absent from connectivity archive '
                f'{connectivity_source}'
            )

        # Richness, HMOD and land use do not depend on connectivity taxon,
        # so bird and mammal members share the same cached inputs.
        richness_key = (sdm, gcm)
        landscape_key = (lum, gcm)

        if richness_key not in richness_cache:
            richness_cache[richness_key] = build_richness_components(
                sdm, gcm, ssp, target,
            )

        if landscape_key not in hmod_cache:
            hmod_cache[landscape_key] = load_hmod(ssp, lum, gcm, target)

        if landscape_key not in landuse_cache:
            landuse_cache[landscape_key] = load_landuse(
                ssp, lum, gcm, target,
            )

        connectivity_class = connectivity_stack[
            connectivity_lookup[identifier]
        ]

        result = common_scale_policy_scores(
            connectivity_class,
            richness_cache[richness_key],
            hmod_cache[landscape_key],
            landuse_cache[landscape_key],
            coverage,
            references,
            target,
        )

        # Store member stacks as float32 to keep memory manageable.
        winners.append(result['winner'].astype(np.float32))
        available_members.append(result['available'])
        connectivity_members.append(connectivity_class)
        alignment_members.append(
            result['recovery_alignment'].astype(np.float32))
        gap_members.append(result['restoration_gap'].astype(np.float32))
        recovery_members.append(result['recovery'].astype(np.float32))

        richness_member = richness_cache[richness_key]
        for name in biodiversity_names:
            if name == 'recipient_opportunity':
                # Use the common availability mask for the action diagnostic.
                values = result[name]
            else:
                values = 100.0 * richness_member[name]
            values_available = np.isfinite(values)
            biodiversity_sums[name][values_available] += values[values_available]
            biodiversity_counts[name][values_available] += 1

        for name in opportunity_names:
            values = result[name]
            values_available = np.isfinite(values)
            opportunity_sums[name][values_available] += values[values_available]
            opportunity_counts[name][values_available] += 1

        absolute_gain = richness_member['gain_absolute_mean']
        absolute_gain_available = np.isfinite(absolute_gain)
        absolute_gain_sum[absolute_gain_available] += absolute_gain[
            absolute_gain_available
        ]
        absolute_gain_count[absolute_gain_available] += 1

        margin = result['winning_margin']
        margin_valid = np.isfinite(margin)
        margin_sum[margin_valid] += margin[margin_valid]
        margin_count[margin_valid] += 1

        for code in POLICY_CODES:
            score_members[code].append(
                result['common_scores'][code].astype(np.float32))
            hotspot = result['hotspot_scores'][code]
            hotspot_valid = np.isfinite(hotspot)
            hotspot_sums[code][hotspot_valid] += hotspot[hotspot_valid]
            hotspot_counts[code][hotspot_valid] += 1

        output_name = (
            f'kmgbf_priority_{ssp}_{taxon}_{sdm}_{gcm}_{lum}_{curve}_'
            f'J1_rad{radius}.tif'
        )

        write_uint8(
            member_output_dir / output_name,
            result['winner'],
            target,
            'Ensemble-member primary KM-GBF policy class',
        )

        inventory.append({
            'member_index': member_index,
            'member_id': identifier,
            'ssp': ssp,
            'radius_pixels': radius,
            'radius_km': radius * 50,
            'connectivity_taxon': taxon,
            'sdm': sdm,
            'gcm': gcm,
            'lum': lum,
            'curve': curve,
            'priority_map': str(member_output_dir / output_name),
        })

        print(
            f'{ssp} rad{radius} [{member_index + 1:03d}/{n_members}] '
            f'{identifier}',
            flush=True,
        )

    winner_stack = np.stack(winners)
    available_stack = np.stack(available_members)
    connectivity_stack_policy = np.stack(connectivity_members)
    alignment_stack = np.stack(alignment_members)
    gap_stack = np.stack(gap_members)
    recovery_stack = np.stack(recovery_members)

    del winners, available_members, connectivity_members
    del alignment_members, gap_members, recovery_members

    score_stacks = {
        code: np.stack(values)
        for code, values in score_members.items()
    }
    del score_members

    available_count = np.sum(available_stack, axis=0)

    def member_mean(stack):
        finite = np.isfinite(stack)
        count = np.sum(finite, axis=0)
        return np.divide(
            np.nansum(stack, axis=0, dtype=np.float64),
            count,
            out=np.full(target['shape'], np.nan, dtype=np.float64),
            where=count > 0,
        )

    # Mean common score for each target across available member scores.
    mean_scores = {
        code: member_mean(stack)
        for code, stack in score_stacks.items()
    }
    mean_hotspots = {
        code: np.divide(
            hotspot_sums[code],
            hotspot_counts[code],
            out=np.full(target['shape'], np.nan, dtype=np.float64),
            where=hotspot_counts[code] > 0,
        )
        for code in POLICY_CODES
    }
    del hotspot_sums, hotspot_counts
    mean_biodiversity = {
        name: np.divide(
            biodiversity_sums[name],
            biodiversity_counts[name],
            out=np.full(target['shape'], np.nan, dtype=np.float64),
            where=biodiversity_counts[name] > 0,
        )
        for name in biodiversity_names
    }
    mean_opportunities = {
        name: np.divide(
            opportunity_sums[name],
            opportunity_counts[name],
            out=np.full(target['shape'], np.nan, dtype=np.float64),
            where=opportunity_counts[name] > 0,
        )
        for name in opportunity_names
    }
    mean_absolute_gain = np.divide(
        absolute_gain_sum,
        absolute_gain_count,
        out=np.full(target['shape'], np.nan, dtype=np.float64),
        where=absolute_gain_count > 0,
    )
    mean_winning_margin = np.divide(
        margin_sum,
        margin_count,
        out=np.full(target['shape'], np.nan, dtype=np.float64),
        where=margin_count > 0,
    )
    del biodiversity_sums, biodiversity_counts
    del opportunity_sums, opportunity_counts
    del absolute_gain_sum, absolute_gain_count
    del margin_sum, margin_count

    # Modal primary class and selection percentages.
    selection_counts = np.stack([
        np.sum((winner_stack == code) & available_stack, axis=0)
        for code in ALL_PRIMARY_CODES
    ])

    max_count = np.max(selection_counts, axis=0)
    tied = selection_counts == max_count[None, :, :]
    n_tied = np.sum(tied, axis=0)

    # Tie-break by mean common-scale score across members.
    # Codes 0 and 97 have no score, so a real policy class wins over them.
    mean_stack = np.stack([
        np.nan_to_num(mean_scores[code], nan=-np.inf)
        if code in mean_scores
        else np.full(target['shape'], -np.inf)
        for code in ALL_PRIMARY_CODES
    ])
    tie_scores = np.where(tied, mean_stack, -np.inf)
    best_tie_score = np.max(tie_scores, axis=0)
    still_tied = np.sum(
        tied & np.isclose(tie_scores, best_tie_score[None, :, :],
                          rtol=0.0, atol=TIE_TOLERANCE),
        axis=0,
    )

    # Without a tie, use the class with the highest count; with a tie,
    # use the tied class with the highest mean score.
    modal_index = np.where(
        n_tied > 1,
        np.argmax(tie_scores, axis=0),
        np.argmax(selection_counts, axis=0),
    )
    modal_class = np.take(
        np.asarray(ALL_PRIMARY_CODES, dtype=np.float64),
        modal_index,
    )

    # Only exact ties in both count and mean score remain unresolved.
    modal_class[(n_tied > 1) & ~np.isfinite(best_tie_score)] = 97
    modal_class[(n_tied > 1) & (still_tied > 1)] = 97
    modal_class[available_count == 0] = np.nan

    modal_agreement = np.divide(
        100.0 * max_count,
        available_count,
        out=np.full(target['shape'], np.nan, dtype=np.float64),
        where=available_count > 0,
    )

    # Connectivity agreement across the same ensemble members,
    # now including connectivity taxon as a source of variation.
    connectivity_valid = np.isin(connectivity_stack_policy, [1, 2, 3, 4])
    connectivity_count = np.sum(connectivity_valid, axis=0)
    connectivity_agreement = {}

    for code, label in CONNECTIVITY_LABELS.items():
        connectivity_agreement[label] = np.divide(
            100.0 * np.sum(
                (connectivity_stack_policy == code) & connectivity_valid,
                axis=0,
            ),
            connectivity_count,
            out=np.full(target['shape'], np.nan, dtype=np.float64),
            where=connectivity_count > 0,
        )

    connectivity_agreement['concentrated'] = (
        connectivity_agreement['channelized']
        + connectivity_agreement['intensified']
    )

    connectivity_total = (
        connectivity_agreement['impeded']
        + connectivity_agreement['diffuse']
        + connectivity_agreement['channelized']
        + connectivity_agreement['intensified']
    )

    if not np.allclose(
        connectivity_total[connectivity_count > 0],
        100.0,
        rtol=0.0,
        atol=1e-5,
    ):
        raise ValueError(
            'Connectivity agreement percentages do not sum to 100')

    mean_alignment = member_mean(alignment_stack)
    mean_gap = member_mean(gap_stack)
    mean_recovery = member_mean(recovery_stack)

    # Modal products.
    write_uint8(
        group_dir / 'kmgbf_modal_priority.tif',
        modal_class,
        target,
        f'Modal KM-GBF policy class across {n_members} ensemble members',
    )

    write_float(
        group_dir / 'kmgbf_modal_priority_agreement_percent.tif',
        modal_agreement,
        target,
        'Percentage of available ensemble members supporting the modal '
        'KM-GBF class',
    )

    write_float(
        group_dir / 'available_member_count.tif',
        available_count.astype(np.float64),
        target,
        'Number of available ensemble members',
    )

    write_float(
        group_dir / 'kmgbf_mean_winning_margin_common_scale.tif',
        mean_winning_margin,
        target,
        'Mean difference between highest and second-highest common target score',
    )

    for name, values in mean_biodiversity.items():
        write_float(
            group_dir / f'{name}_ensemble_mean_0_100.tif',
            values,
            target,
            f'Ensemble mean {name} score on 0-100 scale',
        )

    for name, values in mean_opportunities.items():
        write_float(
            group_dir / f'{name}_ensemble_mean_0_100.tif',
            values,
            target,
            f'Ensemble mean {name} score on 0-100 scale',
        )

    write_float(
        group_dir / 'gain_absolute_mean_richness_ensemble_mean.tif',
        mean_absolute_gain,
        target,
        'Ensemble mean positive richness gain averaged across available taxa',
    )

    # Per-policy selection share and mean continuous score.
    for code in ALL_PRIMARY_CODES:
        frequency = np.divide(
            100.0 * selection_counts[ALL_PRIMARY_CODES.index(code)],
            available_count,
            out=np.full(target['shape'], np.nan, dtype=np.float64),
            where=available_count > 0,
        )

        write_float(
            group_dir / f'policy_selection_percent_code{code}.tif',
            frequency,
            target,
            f'Percentage of available ensemble members selecting policy '
            f'code {code}',
        )

    for code in POLICY_CODES:
        write_float(
            group_dir / f'policy_mean_common_score_code{code}.tif',
            mean_scores[code],
            target,
            f'Mean ensemble-member common-scale score for policy code {code}',
        )
        write_float(
            group_dir / f'policy_mean_hotspot_percentile_code{code}.tif',
            mean_hotspots[code],
            target,
            f'Mean within-target hotspot percentile for policy code {code}',
        )

    for label, values in connectivity_agreement.items():
        write_float(
            group_dir / f'connectivity_{label}_agreement_percent.tif',
            values,
            target,
            f'Percentage of ensemble members classified as {label}',
        )

    write_float(
        group_dir / 'target2_mean_recovery_score.tif',
        mean_recovery,
        target,
        'Mean fixed-reference projected secondary-recovery score',
    )

    write_float(
        group_dir / 'target2_mean_recovery_alignment_index.tif',
        mean_alignment,
        target,
        'Mean Target 2 recovery-alignment index',
    )

    write_float(
        group_dir / 'target2_mean_restoration_gap_index.tif',
        mean_gap,
        target,
        'Mean Target 2 restoration-gap index',
    )

    # Compact member-level score archive for uncertainty attribution.
    archive_path = group_dir / 'kmgbf_member_scores.npz'
    temporary_archive = archive_path.with_name(
        archive_path.stem + '.tmp.npz'
    )

    np.savez_compressed(
        temporary_archive,
        primary_class=np.where(
            np.isfinite(winner_stack),
            winner_stack,
            NODATA_UINT8,
        ).astype(np.uint8),
        available=available_stack.astype(np.uint8),
        connectivity_class=connectivity_stack_policy.astype(np.uint8),
        target2_recovery_alignment=alignment_stack,
        target2_restoration_gap=gap_stack,
        target2_recovery_score=recovery_stack,
        **{
            f'common_score_code_{code}': score_stacks[code]
            for code in POLICY_CODES
        },
        member_id=np.asarray([row['member_id'] for row in inventory]),
        connectivity_taxon=np.asarray(
            [row['connectivity_taxon'] for row in inventory]),
        sdm=np.asarray([row['sdm'] for row in inventory]),
        gcm=np.asarray([row['gcm'] for row in inventory]),
        lum=np.asarray([row['lum'] for row in inventory]),
        curve=np.asarray([row['curve'] for row in inventory]),
        ssp=np.asarray([ssp] * len(inventory)),
        radius=np.asarray([radius] * len(inventory), dtype=np.int16),
    )

    temporary_archive.replace(archive_path)

    inventory_path = group_dir / 'member_inventory.csv'

    with inventory_path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(inventory[0]))
        writer.writeheader()
        writer.writerows(inventory)

    metadata = {
        'ssp': ssp,
        'radius_pixels': radius,
        'radius_km': radius * 50,
        'member_count': n_members,
        'ensemble_design': (
            'connectivity taxon x SDM x GCM x LUM x resistance curve'
        ),
        'direct_path_driven_inputs': True,
        'full_connectivity_archive_member_count': ARCHIVE_MEMBER_COUNT,
        'connectivity_taxa': taxa,
        'members_per_connectivity_taxon': n_members // len(taxa),
        'policy_labels': POLICY_LABELS,
        'connectivity_labels': CONNECTIVITY_LABELS,
        'biodiversity_horizons': {
            'current': (
                'Taxon-balanced present-richness percentile; used for T1, '
                'T8 and T33.'
            ),
            'future': (
                'Taxon-balanced future-richness percentile; used for future '
                'PA/OECM actions T31 and T32.'
            ),
            'forward': (
                'max(current, future); used for T2, T10 and T34 where either '
                'present or future importance justifies action.'
            ),
        },
        'richness_gain_metrics': {
            'source_change': (
                'DeltaR = Rfut - Rpres in richness units, calculated before '
                'any spatial percentile transformation.'
            ),
            'positive_gain': 'Gain = max(DeltaR, 0).',
            'gain_magnitude': (
                'Mean across taxa with positive baseline richness of '
                'clip(Gain/max(Rpres,1),0,1).'
            ),
            'gain_breadth': (
                'Fraction of all paired taxa with Gain > 0; zero-baseline '
                'colonisation can therefore contribute to breadth.'
            ),
            'zero_baseline_rule': (
                'Rpres=0 is excluded from proportional gain magnitude and '
                'decline because a relative rate is undefined, but its '
                'positive absolute gain contributes to gain breadth.'
            ),
            'absolute_gain_diagnostic': (
                'Mean Gain across paired taxa in richness units; not used in '
                'policy scores.'
            ),
            'percentile_difference_used': False,
        },
        'recipient_opportunities': {
            'recipient': 'GM(B_future, gain magnitude, gain breadth)',
            'expansion': 'GM(recipient, 1-PA/OECM coverage)',
            'restoration': (
                'GM(recipient, HMOD), available only for impeded or diffuse '
                'connectivity'
            ),
            'role_in_winner_selection': False,
        },
        'climate_vulnerability_formula': (
            'M = absolute mean proportional richness decline across available '
            'taxa; Q = fraction of available taxa with positive decline; '
            'vulnerability = sqrt(M * Q). M is not spatially reranked.'
        ),
        'scaling_reference_file': str(SCALING_REFERENCE_PATH),
        'scaling_references_p99': references,
        'continuous_inputs': {
            'biodiversity_current': '0-1 taxon-balanced spatial percentile',
            'biodiversity_future': '0-1 taxon-balanced spatial percentile',
            'biodiversity_forward': 'max(current, future)',
            'gain_magnitude': '0-1 capped mean proportional richness gain',
            'gain_breadth': '0-1 fraction of taxa with positive richness gain',
            'recipient_opportunity': (
                '0-1 GM(future biodiversity, gain magnitude, gain breadth)'
            ),
            'climate_vulnerability': '0-1 absolute decline/breadth index',
            'HMOD': 'native 0-1 value',
            'PA_OECM_coverage': 'native terrestrial fraction 0-1',
            'coverage_gap': '1 - PA/OECM fraction',
            'primary_loss': 'raw positive change / pooled P99, clipped 0-1',
            'conversion': 'raw positive change / pooled P99, clipped 0-1',
            'production_expansion': (
                'raw production_magnitude / pooled P99, clipped 0-1'
            ),
            'future_production_extent': 'native 0-1 cropland+pasture fraction',
            'secondary_recovery': (
                'raw positive change / pooled P99, clipped 0-1'
            ),
        },
        'categorical_input': (
            'Only connectivity state is categorical: impeded, diffuse, '
            'channelized, intensified.'
        ),
        'target_scores': {
            '1': 'GM(B_current, primary_loss, vulnerability)',
            '21': (
                'GM(B_forward, HMOD); available only where connectivity is '
                'impeded or diffuse'
            ),
            '8': 'GM(B_current, vulnerability)',
            '10': (
                'GM(B_forward, future production extent, production expansion*)'
            ),
            '31': (
                'GM(B_future, 1-HMOD, 1-PA/OECM, persistence); available '
                'only where connectivity is diffuse'
            ),
            '32': (
                'GM(B_future, 1-PA/OECM); available only where connectivity '
                'is channelized or intensified'
            ),
            '33': 'GM(B_current, PA/OECM coverage, vulnerability)',
            '34': 'GM(B_forward, PA/OECM coverage, conversion*)',
        },
        'winner_selection': (
            'Winner is the maximum 100 * raw geometric-mean score. Final '
            'target scores are not spatially reranked.'
        ),
        'hotspot_diagnostics': (
            'Positive raw scores are spatially percentile-ranked separately '
            'within each target only for diagnostic hotspot maps; these '
            'percentiles never determine winners.'
        ),
        'biodiversity_diagnostics': (
            'Ensemble-mean current, future, forward, gain magnitude, gain '
            'breadth and recipient-opportunity maps are written on a 0-100 '
            'scale. Absolute positive gain is written separately in richness '
            'units. Gain metrics come from Rpres and Rfut, not from differences '
            'between spatial percentiles.'
        ),
        'availability_rule': (
            'All eight competing target scores use one shared cell-by-member '
            'availability mask. Recipient diagnostics additionally require a '
            'finite RecipientOpportunity; missing gain data do not suppress '
            'otherwise available policy scores.'
        ),
        'winning_margin_diagnostic': (
            'Mean common-scale difference between the highest and second-highest '
            'target score; it does not affect winner selection.'
        ),
        'modal_tie_break': (
            'Classes tied for the highest member count are resolved by the '
            'highest mean common-scale score; exact remaining ties = 97. '
            'Cells where all available targets score zero are coded 0.'
        ),
        'target2_recovery_context': {
            'alignment': 'T2 common score * fixed-reference recovery score',
            'gap': 'T2 common score * (1 - fixed-reference recovery score)',
            'recovery_never_suppresses_need': True,
        },
        'pa_oecm_raster': str(PA_OECM_PATH),
        'connectivity_class_archive': str(connectivity_source),
        'member_score_archive': str(archive_path),
        'member_inventory': str(inventory_path),
        'interpretation': (
            'Selection percentages and connectivity agreement describe the '
            'specified model ensemble and are not calibrated probabilities.'
        ),
    }

    (group_dir / 'metadata.json').write_text(
        json.dumps(metadata, indent=2),
        encoding='utf-8',
    )

    print(f'COMPLETED {ssp} | radius {radius} | {n_members} members',
          flush=True)

# =============================================================================
# MAIN
# =============================================================================


def main():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('--ssp', nargs='+', choices=SSPS, default=SSPS)
    parser.add_argument(
        '--radius', nargs='+', type=int, choices=RADII, default=RADII,
    )
    parser.add_argument(
        '--connectivity-taxa',
        nargs='+',
        choices=CONNECTIVITY_TAXA,
        default=CONNECTIVITY_TAXA,
        help='Default uses both taxa (180 members). A single taxon gives '
             'a 90-member sensitivity run in a separate output folder.',
    )

    args = parser.parse_args()

    taxa = [t for t in CONNECTIVITY_TAXA if t in args.connectivity_taxa]

    if taxa == CONNECTIVITY_TAXA:
        output_root = OUTPUT_ROOT_ALL_TAXA
    else:
        output_root = BASE / (
            f'kmgbf_memberwise_common_scale_{"_".join(taxa)}_only_'
            'recipient_aware'
        )

    references, _ = load_scaling_references()
    target = load_template()
    coverage = read_raster(
        PA_OECM_PATH,
        target,
        minimum=0.0,
        maximum=1.0,
    )

    output_root.mkdir(parents=True, exist_ok=True)
    print(f'Connectivity taxa: {taxa} | output: {output_root}', flush=True)
    print(f'Common P99 references: {references}', flush=True)

    for ssp in args.ssp:
        for radius in args.radius:
            process_group(
                ssp, radius, target, coverage, references, taxa, output_root,
            )

    print('ALL REQUESTED COMMON-SCALE MEMBER-WISE PRIORITISATIONS COMPLETED')


if __name__ == '__main__':
    main()
