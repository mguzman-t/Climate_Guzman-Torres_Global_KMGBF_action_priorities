#!/usr/bin/env python3
"""
Calibrate KM-GBF target scores using fixed empirical percentile distributions.

The script reads the existing "kmgbf_member_scores.npz" archives and
constructs a separate calibration distribution for each strategic target.
For a given target, all finite, positive and eligible scores are pooled across
the requested SSPs, connectivity radii, ensemble members and grid cells.

Each target therefore has its own pooled distribution of positive scores.
These distributions are constructed once and then held fixed. The same
target-specific calibration is applied to every SSP and connectivity radius,
so calibrated scores are directly comparable across scenario-radius
combinations.

For each positive score, calibration reports its empirical percentile within
the pooled distribution of positive scores for the same target. For example,
a calibrated value of 80 indicates that the score is higher than approximately
80% of the positive scores in that target's pooled calibration distribution.
The calibrated value is therefore a relative rank within a target.


Important considerations
------------------------

* The original common-range scores are retained unchanged. Calibration produces
  additional scores and does not overwrite the underlying ecological values.

* Only finite and eligible observations contribute to calibration. For targets
  with a connectivity-class eligibility rule, cells outside the applicable
  connectivity classes are excluded from the calibration distribution rather
  than represented by zeros.

* Eligible scores equal to zero remain zero after calibration. They are retained
  in audit counts but are excluded from the positive-score distribution used
  to calculate empirical percentiles. Consequently, a member for which all
  eligible target scores are zero remains classified as having no positive
  action.

* The empirical cumulative distribution function (ECDF) is calculated
  separately for each target using its pooled positive scores. Positive scores
  are converted to empirical mid-rank percentiles, with tied raw values
  receiving the same mid-rank percentile.

* Member-level target selection is recalculated using the calibrated scores.
  Consequently, winning actions, winning margins, member priority maps,
  action-selection percentages, modal priority and modal agreement are all
  recalculated after calibration.

* Within-group hotspot percentiles are retained only as diagnostic outputs.
  They are not used to select the winning strategic action.

* Target 2 retains internal policy code 21. Unresolved positive ties retain
  internal code 97.

* Recipient-opportunity, recovery and other diagnostics that do not depend on
  winner selection are retained unchanged where they are available in the
  source outputs. These diagnostics do not determine the winning action.


Calibration implementation
--------------------------

The calibration uses the observed scores directly rather than approximating
their distributions with histograms. For each target, all positive eligible
scores are pooled and sorted into a memory-mapped ".npy" array.

A positive score is then located within this sorted distribution and converted
to its empirical mid-rank percentile. The resulting sorted arrays constitute
the fixed calibration references. They can be inspected and reused in later
runs without reconstructing the pooled distributions.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import rasterio

# -----------------------------------------------------------------------------
# Fixed scientific identifiers
# -----------------------------------------------------------------------------

POLICY_CODES = (1, 21, 8, 10, 31, 32, 33, 34)
ALL_PRIMARY_CODES = (0, 1, 21, 8, 10, 31, 32, 33, 34, 97)
POLICY_LABELS = {
    0: 'No positive action',
    1: 'T1 spatial planning',
    21: 'T2 restoration',
    8: 'T8 climate adaptation',
    10: 'T10 production landscapes',
    31: 'T31 maintain diffuse landscapes',
    32: 'T32 safeguard concentrated connectivity',
    33: 'T33 support PA/OECM under climate vulnerability',
    34: 'T34 safeguard PA/OECM under land-use pressure',
    97: 'Unresolved positive tie',
}

SSPS = ('ssp126', 'ssp370', 'ssp585')
RADII = (6, 14, 20)
SOURCE_ROOT = Path(
    '/capacity/occr_davin/mguzman/P1/Connectivity/approach2/'
    'kmgbf_memberwise_common_scale_180'
)
OUTPUT_ROOT = Path(
    '/capacity/occr_davin/mguzman/P1/Connectivity/approach2/'
    'kmgbf_memberwise_t_calibration_positive_180'
)

NODATA_UINT8 = np.uint8(255)
NODATA_FLOAT32 = np.float32(-9999.0)
RAW_SCORE_MIN = 0.0
RAW_SCORE_MAX = 100.0
RAW_RANGE_TOLERANCE = 1.0e-3
TIE_TOLERANCE = 1.0e-6
CHUNK_CELLS = 1_000_000
CALIBRATION_VERSION = 2
CALIBRATION_METHOD = (
    'exact positive-only pooled target-specific empirical mid-rank ECDF'
)
CALIBRATION_ZERO_RULE = (
    'Eligible raw scores equal to zero map to calibrated zero and are retained '
    'in audit counts, but are excluded from the target-specific ECDF reference '
    'distribution. Positive scores are calibrated only against positive '
    'eligible pooled scores.'
)

# Archive arrays copied to the calibrated archive without changing their
# scientific meaning.  The first three preserve the requested recovery
# diagnostics.  Additional member descriptors are handled separately.
PASSTHROUGH_ARRAY_KEYS = (
    'target2_recovery_alignment',
    'target2_restoration_gap',
    'target2_recovery_score',
)
MEMBER_DESCRIPTOR_KEYS = (
    'member_id',
    'connectivity_taxon',
    'sdm',
    'gcm',
    'lum',
    'curve',
    'ssp',
    'radius',
)

# These products depend on the old, uncalibrated winners and must not be copied
# into the new output group.  They are recalculated below.
RECALCULATED_NAMES = {
    'kmgbf_modal_priority.tif',
    'kmgbf_modal_priority_agreement_percent.tif',
    'available_member_count.tif',
    'kmgbf_mean_winning_margin_common_scale.tif',
    'kmgbf_mean_winning_margin_fixed_pooled_calibrated.tif',
    'kmgbf_member_scores_fixed_pooled_calibrated.npz',
    'member_inventory.csv',
    'metadata.json',
    'raw_archive_source.txt',
}


@dataclass(frozen=True)
class Group:
    """A single SSP-by-radius source group."""

    ssp: str
    radius: int
    directory: Path
    archive: Path

    @property
    def name(self) -> str:
        return f'{self.ssp}_rad{self.radius}'


def utc_now() -> str:
    """Return a reproducible ISO-style UTC timestamp for provenance."""
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def atomic_write_json(path: Path, value: Mapping) -> None:
    """Write JSON through a temporary file, then replace atomically."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True),
        encoding='utf-8',
    )
    temporary.replace(path)


def sha256_file(path: Path, chunk_bytes: int = 16 * 1024 * 1024) -> str:
    """Calculate a file digest without reading the whole file into memory."""
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        while True:
            block = handle.read(chunk_bytes)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def discover_groups(
    root: Path,
    ssps: Sequence[str],
    radii: Sequence[int],
    expected_members: int | None,
) -> list[Group]:
    """Resolve and validate all requested source score archives."""
    groups: list[Group] = []
    missing: list[Path] = []

    for ssp in ssps:
        for radius in radii:
            directory = root / f'{ssp}_rad{radius}'
            archive = directory / 'kmgbf_member_scores.npz'
            if not archive.is_file():
                missing.append(archive)
                continue
            groups.append(Group(ssp, radius, directory, archive))

    if missing:
        formatted = '\n'.join(f'  - {path}' for path in missing)
        raise FileNotFoundError(
            'The fixed calibration requires every requested source archive. '
            f'Missing:\n{formatted}'
        )

    for group in groups:
        with np.load(group.archive, allow_pickle=False) as source:
            validate_archive_structure(source, group, expected_members)

    return groups


def validate_archive_structure(
    source: Mapping[str, np.ndarray],
    group: Group,
    expected_members: int | None,
) -> tuple[int, int, int]:
    """Check archive keys and the common member-by-row-by-column shape."""
    required = {
        'available',
        'connectivity_class',
        *{f'common_score_code_{code}' for code in POLICY_CODES},
    }
    missing = sorted(required.difference(source.keys()))
    if missing:
        raise KeyError(f'{group.archive} is missing arrays: {missing}')

    shape = tuple(source['available'].shape)
    if len(shape) != 3:
        raise ValueError(
            f'{group.archive}: available must be [member, row, column], '
            f'not {shape}'
        )
    if expected_members is not None and shape[0] != expected_members:
        raise ValueError(
            f'{group.archive}: expected {expected_members} members, '
            f'found {shape[0]}'
        )

    for key in required.difference({'available'}):
        if tuple(source[key].shape) != shape:
            raise ValueError(
                f'{group.archive}: {key} has shape {source[key].shape}; '
                f'expected {shape}'
            )

    return shape


def connectivity_gate(code: int, connectivity: np.ndarray) -> np.ndarray:
    """Return the member-specific categorical eligibility gate for a target."""
    if code == 21:  # T2 restoration: impeded or diffuse
        return np.isin(connectivity, (1, 2))
    if code == 31:  # T31: diffuse
        return connectivity == 2
    if code == 32:  # T32: channelised or intensified
        return np.isin(connectivity, (3, 4))
    return np.ones(connectivity.shape, dtype=bool)


def eligible_score_mask(
    code: int,
    raw_score: np.ndarray,
    available: np.ndarray,
    connectivity: np.ndarray,
) -> np.ndarray:
    """
    Identify finite observations eligible for a target.

    Applying the connectivity gate here is intentional even when old archives
    already contain NaN outside the gate.  It prevents an archive containing
    zero-filled ineligible cells from contaminating the pooled distribution.
    """
    return (
        available
        & connectivity_gate(code, connectivity)
        & np.isfinite(raw_score)
    )


def checked_raw_scores(
    values: np.ndarray,
    mask: np.ndarray,
    source_description: str,
) -> np.ndarray:
    """Validate the raw 0-100 range and clip only floating-point overshoot."""
    selected = np.asarray(values[mask], dtype=np.float32)
    if selected.size == 0:
        return selected

    minimum = float(np.min(selected))
    maximum = float(np.max(selected))
    if minimum < RAW_SCORE_MIN - RAW_RANGE_TOLERANCE:
        raise ValueError(
            f'{source_description}: minimum eligible score {minimum} is below 0'
        )
    if maximum > RAW_SCORE_MAX + RAW_RANGE_TOLERANCE:
        raise ValueError(
            f'{source_description}: maximum eligible score {maximum} exceeds 100'
        )

    return np.clip(selected, RAW_SCORE_MIN, RAW_SCORE_MAX)


def sorted_quantile(sorted_values: np.ndarray, probability: float) -> float:
    """Linear quantile from an already sorted one-dimensional array."""
    if sorted_values.size == 0:
        return float('nan')
    position = probability * (sorted_values.size - 1)
    lower = int(np.floor(position))
    upper = int(np.ceil(position))
    if lower == upper:
        return float(sorted_values[lower])
    fraction = position - lower
    return float(
        sorted_values[lower] * (1.0 - fraction)
        + sorted_values[upper] * fraction
    )


def build_fixed_calibrations(
    groups: Sequence[Group],
    calibration_dir: Path,
    overwrite: bool,
) -> dict:
    """
    Build exact positive-only pooled ECDF references, one target at a time.

    Memory-mapped arrays make the peak memory independent of the total number
    of pooled observations.  Disk use is approximately four bytes per positive
    eligible observation per target because references are stored as sorted
    float32.
    """
    calibration_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = calibration_dir / 'calibration_manifest.json'

    if manifest_path.exists() and not overwrite:
        raise FileExistsError(
            f'{manifest_path} already exists. Use --reuse-calibration to apply '
            'it, or --overwrite-calibration to rebuild it explicitly.'
        )

    group_records = [
        {
            'ssp': group.ssp,
            'radius_pixels': group.radius,
            'source_archive': str(group.archive.resolve()),
            'source_archive_sha256': sha256_file(group.archive),
        }
        for group in groups
    ]

    targets: dict[str, dict] = {}

    for code in POLICY_CODES:
        print(
            f'CALIBRATION T{2 if code == 21 else code}: counting', flush=True)
        eligible_count = 0
        positive_count = 0
        zero_count = 0

        for group in groups:
            with np.load(group.archive, allow_pickle=False) as source:
                available = np.asarray(source['available'], dtype=bool)
                connectivity = np.asarray(source['connectivity_class'])
                raw = np.asarray(
                    source[f'common_score_code_{code}'], dtype=np.float32,
                )
                mask = eligible_score_mask(
                    code, raw, available, connectivity,
                )
                selected = checked_raw_scores(
                    raw, mask, f'{group.name}, policy code {code}',
                )
                positive = selected[selected > 0.0]
                eligible_count += int(selected.size)
                positive_count += int(positive.size)
                zero_count += int(np.count_nonzero(selected <= 0.0))

        if eligible_count == 0:
            raise ValueError(
                f'Policy code {code} has no finite eligible observations '
                'across the calibration groups'
            )
        if positive_count == 0:
            raise ValueError(
                f'Policy code {code} has no positive eligible observations '
                'across the calibration groups'
            )

        reference_path = calibration_dir / (
            f'policy_code_{code}_pooled_sorted_positive_raw_scores.npy'
        )
        if reference_path.exists():
            if not overwrite:
                raise FileExistsError(reference_path)
            reference_path.unlink()

        print(
            f'CALIBRATION code {code}: {eligible_count:,} eligible; '
            f'{positive_count:,} positive reference values; '
            f'{zero_count:,} audited zeros',
            flush=True,
        )
        reference = np.lib.format.open_memmap(
            reference_path,
            mode='w+',
            dtype=np.float32,
            shape=(positive_count,),
        )
        cursor = 0

        for group in groups:
            with np.load(group.archive, allow_pickle=False) as source:
                available = np.asarray(source['available'], dtype=bool)
                connectivity = np.asarray(source['connectivity_class'])
                raw = np.asarray(
                    source[f'common_score_code_{code}'], dtype=np.float32,
                )
                mask = eligible_score_mask(
                    code, raw, available, connectivity,
                )
                selected = checked_raw_scores(
                    raw, mask, f'{group.name}, policy code {code}',
                )
                positive = selected[selected > 0.0]
                next_cursor = cursor + positive.size
                reference[cursor:next_cursor] = positive
                cursor = next_cursor

        if cursor != positive_count:
            raise RuntimeError(
                f'Policy code {code}: counted {positive_count} positive '
                f'observations but wrote {cursor}'
            )

        print(f'CALIBRATION code {code}: sorting', flush=True)
        reference.sort(kind='quicksort')
        reference.flush()

        targets[str(code)] = {
            'policy_label': POLICY_LABELS[code],
            'reference_file': reference_path.name,
            'eligible_observation_count': eligible_count,
            'eligible_zero_count': zero_count,
            'eligible_zero_percent': (
                100.0 * zero_count / eligible_count
            ),
            'positive_reference_count': positive_count,
            'positive_reference_percent_of_eligible': (
                100.0 * positive_count / eligible_count
            ),
            'reference_excludes_zeros': True,
            'minimum_positive': float(reference[0]),
            'p01_positive': sorted_quantile(reference, 0.01),
            'p05_positive': sorted_quantile(reference, 0.05),
            'p25_positive': sorted_quantile(reference, 0.25),
            'median_positive': sorted_quantile(reference, 0.50),
            'p75_positive': sorted_quantile(reference, 0.75),
            'p95_positive': sorted_quantile(reference, 0.95),
            'p99_positive': sorted_quantile(reference, 0.99),
            'maximum_positive': float(reference[-1]),
        }
        del reference

    manifest = {
        'calibration_version': CALIBRATION_VERSION,
        'created_utc': utc_now(),
        'method': CALIBRATION_METHOD,
        'raw_score_scale': [RAW_SCORE_MIN, RAW_SCORE_MAX],
        'calibrated_score_scale': [0.0, 100.0],
        'reference_population': (
            'Finite raw scores greater than zero that satisfy available=1 and '
            'the target-specific connectivity eligibility gate.'
        ),
        'midrank_formula': (
            '100 * (count(reference < x) + 0.5 * count(reference == x)) / '
            'positive_reference_count, for eligible x > 0'
        ),
        'zero_rule': CALIBRATION_ZERO_RULE,
        'decision_interpretation': (
            'For targets with positive evidence, the calibrated score '
            'represents the empirical standing of that positive evidence '
            'within the fixed positive-score distribution for the same target, '
            'pooled across all requested SSPs, radii and ensemble members.'
        ),
        'eligibility_rule': (
            'Only finite scores with available=1 and the target-specific '
            'connectivity gate satisfied are pooled. Connectivity-ineligible '
            'cells are excluded, not represented as zero.'
        ),
        'connectivity_gates': {
            '1': 'no additional gate',
            '21': 'connectivity class 1 or 2',
            '8': 'no additional gate',
            '10': 'no additional gate',
            '31': 'connectivity class 2',
            '32': 'connectivity class 3 or 4',
            '33': 'no additional gate',
            '34': 'no additional gate',
        },
        'frozen_across_ssps_and_radii': True,
        'calibration_groups': group_records,
        'targets': targets,
    }
    atomic_write_json(manifest_path, manifest)
    print(f'WROTE {manifest_path}', flush=True)
    return manifest


def load_fixed_calibrations(calibration_dir: Path) -> tuple[dict, dict[int, np.ndarray]]:
    """Open frozen sorted references as read-only memory maps."""
    manifest_path = calibration_dir / 'calibration_manifest.json'
    if not manifest_path.is_file():
        raise FileNotFoundError(manifest_path)
    manifest = json.loads(manifest_path.read_text(encoding='utf-8'))

    if manifest.get('calibration_version') != CALIBRATION_VERSION:
        raise ValueError(
            f'Unsupported calibration version in {manifest_path}: '
            f"{manifest.get('calibration_version')}"
        )
    if manifest.get('method') != CALIBRATION_METHOD:
        raise ValueError(
            f'Unsupported calibration method in {manifest_path}: '
            f"{manifest.get('method')!r}"
        )

    references: dict[int, np.ndarray] = {}
    for code in POLICY_CODES:
        record = manifest.get('targets', {}).get(str(code))
        if not record:
            raise KeyError(
                f'{manifest_path}: no target record for code {code}')
        path = calibration_dir / record['reference_file']
        reference = np.load(path, mmap_mode='r', allow_pickle=False)
        if reference.ndim != 1 or reference.size == 0:
            raise ValueError(
                f'{path}: expected a non-empty one-dimensional array')
        if reference.size != record['positive_reference_count']:
            raise ValueError(
                f'{path}: manifest positive-reference count is '
                f"{record['positive_reference_count']}, file count is "
                f'{reference.size}'
            )
        if not record.get('reference_excludes_zeros', False):
            raise ValueError(
                f'{manifest_path}: policy code {code} is not marked as using '
                'a zero-excluding reference'
            )
        if record['eligible_observation_count'] != (
            record['eligible_zero_count'] + record['positive_reference_count']
        ):
            raise ValueError(
                f'{manifest_path}: policy code {code} audit counts are '
                'internally inconsistent'
            )
        previous = -np.inf
        for start in range(0, reference.size, CHUNK_CELLS):
            stop = min(start + CHUNK_CELLS, reference.size)
            block = np.asarray(reference[start:stop])
            if not np.all(np.isfinite(block)):
                raise ValueError(
                    f'{path}: reference contains non-finite values')
            if np.any(block <= 0.0):
                raise ValueError(
                    f'{path}: positive-only reference contains a value less '
                    'than or equal to zero'
                )
            if float(block[0]) < previous or np.any(block[1:] < block[:-1]):
                raise ValueError(f'{path}: reference is not sorted ascending')
            previous = float(block[-1])
        references[code] = reference

    return manifest, references


def calibrate_scores(
    raw_score: np.ndarray,
    eligible: np.ndarray,
    sorted_reference: np.ndarray,
    chunk_cells: int,
) -> np.ndarray:
    """
    Map eligible raw scores to fixed positive-only empirical mid-ranks.

    Chunking bounds the temporary memory used by two ``searchsorted`` index
    arrays.  Ineligible observations remain NaN.  Eligible zeros remain zero
    and never enter the positive reference distribution.
    """
    if raw_score.shape != eligible.shape:
        raise ValueError('raw_score and eligible must have the same shape')

    output = np.full(raw_score.shape, np.nan, dtype=np.float32)
    raw_flat = raw_score.reshape(-1)
    eligible_flat = eligible.reshape(-1)
    output_flat = output.reshape(-1)
    n_reference = sorted_reference.size
    if n_reference == 0:
        raise ValueError('Positive-only calibration reference is empty')

    for start in range(0, raw_flat.size, chunk_cells):
        stop = min(start + chunk_cells, raw_flat.size)
        local_eligible = eligible_flat[start:stop]
        if not np.any(local_eligible):
            continue

        local_raw = np.asarray(raw_flat[start:stop], dtype=np.float32)
        positions = np.flatnonzero(local_eligible)
        values = np.clip(
            local_raw[positions], RAW_SCORE_MIN, RAW_SCORE_MAX,
        )

        zero = values <= 0.0
        calibrated = np.zeros(values.shape, dtype=np.float32)
        if np.any(~zero):
            positive = values[~zero]
            left = np.searchsorted(sorted_reference, positive, side='left')
            right = np.searchsorted(sorted_reference, positive, side='right')
            # left = number strictly below; right - left = number tied.
            # Therefore 0.5 * (left + right) is the exact empirical mid-rank.
            calibrated[~zero] = (
                50.0 * (left.astype(np.float64) + right.astype(np.float64))
                / n_reference
            ).astype(np.float32)

        output_flat[start + positions] = calibrated

    return output


def member_mean(stack: np.ndarray) -> np.ndarray:
    """Mean over finite ensemble-member values, retaining NaN if none exist."""
    finite = np.isfinite(stack)
    count = np.sum(finite, axis=0)
    return np.divide(
        np.nansum(stack, axis=0, dtype=np.float64),
        count,
        out=np.full(stack.shape[1:], np.nan, dtype=np.float64),
        where=count > 0,
    )


def find_template_raster(group: Group) -> Path:
    """Use an existing group output solely for georeferencing and dimensions."""
    preferred = group.directory / 'kmgbf_modal_priority.tif'
    if preferred.is_file():
        return preferred
    candidates = sorted(group.directory.glob('*.tif'))
    if not candidates:
        raise FileNotFoundError(
            f'{group.directory}: no GeoTIFF is available as an output template'
        )
    return candidates[0]


def template_profile(path: Path, expected_shape: tuple[int, int]) -> dict:
    """Read and validate a raster profile without using its scientific values."""
    with rasterio.open(path) as source:
        if (source.height, source.width) != expected_shape:
            raise ValueError(
                f'{path}: raster shape {(source.height, source.width)} does not '
                f'match archive shape {expected_shape}'
            )
        return source.profile.copy()


def write_uint8(path: Path, values: np.ndarray, profile: Mapping, description: str) -> None:
    """Write a categorical map with 255 as nodata."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if values.dtype != np.uint8:
        output = np.where(np.isfinite(values), values,
                          NODATA_UINT8).astype(np.uint8)
    else:
        output = values

    raster_profile = dict(profile)
    # Very small synthetic or regional rasters cannot use TIFF tile dimensions
    # that are multiples of 16.  Strip inherited block sizes and use strips for
    # those cases; global production rasters remain tiled.
    raster_profile.pop('blockxsize', None)
    raster_profile.pop('blockysize', None)
    use_tiles = values.shape[-2] >= 16 and values.shape[-1] >= 16
    raster_profile.update(
        driver='GTiff',
        count=1,
        dtype='uint8',
        nodata=int(NODATA_UINT8),
        compress='deflate',
        predictor=2,
        tiled=use_tiles,
        BIGTIFF='IF_SAFER',
    )
    with rasterio.open(path, 'w', **raster_profile) as destination:
        destination.write(output, 1)
        destination.set_band_description(1, description)
        destination.update_tags(
            calibration_version=str(CALIBRATION_VERSION),
            calibration=CALIBRATION_METHOD,
            score_role=(
                'winner derived from fixed positive-only pooled '
                'target-specific calibration'
            ),
            policy_code_21='T2',
            unresolved_positive_tie_code='97',
        )


def write_float(path: Path, values: np.ndarray, profile: Mapping, description: str) -> None:
    """Write a continuous float32 diagnostic with explicit nodata."""
    path.parent.mkdir(parents=True, exist_ok=True)
    output = np.where(
        np.isfinite(values), values, NODATA_FLOAT32,
    ).astype(np.float32)
    raster_profile = dict(profile)
    raster_profile.pop('blockxsize', None)
    raster_profile.pop('blockysize', None)
    use_tiles = values.shape[-2] >= 16 and values.shape[-1] >= 16
    raster_profile.update(
        driver='GTiff',
        count=1,
        dtype='float32',
        nodata=float(NODATA_FLOAT32),
        compress='deflate',
        predictor=3,
        tiled=use_tiles,
        BIGTIFF='IF_SAFER',
    )
    with rasterio.open(path, 'w', **raster_profile) as destination:
        destination.write(output, 1)
        destination.set_band_description(1, description)
        destination.update_tags(
            calibration_version=str(CALIBRATION_VERSION),
            calibration=CALIBRATION_METHOD,
            calibration_zero_treatment=(
                'eligible raw zeros remain zero and are excluded from ECDF '
                'references'
            ),
        )


def should_recalculate(name: str) -> bool:
    """Return True for an old winner-dependent product."""
    return (
        name in RECALCULATED_NAMES
        or name.startswith('policy_selection_percent_code')
        or name.startswith('member_priority_maps')
    )


def preserve_existing_products(
    source_dir: Path,
    destination_dir: Path,
    mode: str,
) -> list[str]:
    """
    Retain source diagnostics that do not depend on winner selection.

    Hotspot, recipient-opportunity, recovery, biodiversity and connectivity
    diagnostics are therefore retained.  Old winner products are excluded and
    recalculated.  The large raw archive is not duplicated; its absolute path
    is recorded separately and its raw means are written explicitly.
    """
    if mode == 'none':
        return []

    preserved: list[str] = []
    for source in sorted(source_dir.iterdir()):
        name = source.name
        if name == 'kmgbf_member_scores.npz' or should_recalculate(name):
            continue
        destination = destination_dir / name
        if destination.exists() or destination.is_symlink():
            continue

        if source.is_dir():
            # No source directory other than member_priority_maps is currently
            # expected, but preserve future diagnostic directories sensibly.
            if mode == 'copy':
                shutil.copytree(source, destination)
            elif mode == 'symlink':
                destination.symlink_to(
                    source.resolve(), target_is_directory=True)
            elif mode == 'hardlink':
                shutil.copytree(source, destination, copy_function=os.link)
        else:
            destination.parent.mkdir(parents=True, exist_ok=True)
            if mode == 'copy':
                shutil.copy2(source, destination)
            elif mode == 'symlink':
                destination.symlink_to(source.resolve())
            elif mode == 'hardlink':
                try:
                    os.link(source, destination)
                except OSError:
                    shutil.copy2(source, destination)
        preserved.append(name)

    return preserved


def load_inventory(
    group: Group,
    source: Mapping[str, np.ndarray],
    n_members: int,
) -> list[dict[str, str]]:
    """Read the source inventory or reconstruct it from archive descriptors."""
    path = group.directory / 'member_inventory.csv'
    if path.is_file():
        with path.open(newline='', encoding='utf-8') as handle:
            rows = list(csv.DictReader(handle))
        if len(rows) != n_members:
            raise ValueError(
                f'{path}: expected {n_members} rows, found {len(rows)}'
            )
        return rows

    def strings(key: str, fallback: str) -> list[str]:
        if key not in source:
            return [fallback] * n_members
        values = np.asarray(source[key]).astype(str).tolist()
        if len(values) != n_members:
            raise ValueError(
                f'{group.archive}: {key} has {len(values)} entries')
        return values

    member_ids = strings('member_id', 'unknown_member')
    taxa = strings('connectivity_taxon', 'unknown_taxon')
    sdms = strings('sdm', 'unknown_sdm')
    gcms = strings('gcm', 'unknown_gcm')
    lums = strings('lum', 'unknown_lum')
    curves = strings('curve', 'unknown_curve')

    rows: list[dict[str, str]] = []
    for index in range(n_members):
        filename = (
            f'kmgbf_priority_{group.ssp}_{taxa[index]}_{sdms[index]}_'
            f'{gcms[index]}_{lums[index]}_{curves[index]}_J1_'
            f'rad{group.radius}.tif'
        )
        rows.append({
            'member_index': str(index),
            'member_id': member_ids[index],
            'ssp': group.ssp,
            'radius_pixels': str(group.radius),
            'radius_km': str(group.radius * 50),
            'connectivity_taxon': taxa[index],
            'sdm': sdms[index],
            'gcm': gcms[index],
            'lum': lums[index],
            'curve': curves[index],
            'priority_map': filename,
        })
    return rows


def member_map_name(row: Mapping[str, str], index: int) -> str:
    """Retain the old map basename where available."""
    old_path = row.get('priority_map', '')
    if old_path:
        return Path(old_path).name
    member_id = row.get('member_id', f'member_{index:03d}')
    return f'kmgbf_priority_{member_id}.tif'


def save_inventory(path: Path, rows: Sequence[Mapping[str, str]]) -> None:
    """Write a member inventory with paths to the recalibrated maps."""
    if not rows:
        raise ValueError('Cannot write an empty inventory')
    fieldnames = list(rows[0].keys())
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def update_top_two(
    best: np.ndarray,
    second: np.ndarray,
    best_code: np.ndarray,
    candidate: np.ndarray,
    code: int,
) -> None:
    """Update the highest and second-highest finite calibrated scores in place."""
    filled = np.where(np.isfinite(candidate), candidate, -np.inf)
    previous_best = best.copy()
    better = filled > previous_best
    second[...] = np.where(
        better,
        previous_best,
        np.maximum(second, filled),
    )
    best[...] = np.where(better, filled, previous_best)
    best_code[better] = code


def compute_modal_products(
    winner: np.ndarray,
    member_available: np.ndarray,
    mean_calibrated: Mapping[int, np.ndarray],
    tie_tolerance: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Calculate selection counts, modal class and modal agreement."""
    shape = winner.shape[1:]
    available_count = np.sum(member_available, axis=0, dtype=np.uint16)
    selection_counts = np.stack([
        np.sum((winner == code) & member_available, axis=0, dtype=np.uint16)
        for code in ALL_PRIMARY_CODES
    ])
    max_count = np.max(selection_counts, axis=0)
    tied = selection_counts == max_count[None, :, :]
    number_tied = np.sum(tied, axis=0)

    mean_stack = np.stack([
        np.nan_to_num(mean_calibrated[code], nan=-np.inf)
        if code in mean_calibrated
        else np.full(shape, -np.inf, dtype=np.float64)
        for code in ALL_PRIMARY_CODES
    ])
    tie_scores = np.where(tied, mean_stack, -np.inf)
    best_tie_score = np.max(tie_scores, axis=0)
    best_tie_index = np.argmax(tie_scores, axis=0)
    score_tie_count = np.sum(
        tied
        & np.isfinite(tie_scores)
        & np.isclose(
            tie_scores,
            best_tie_score[None, :, :],
            rtol=0.0,
            atol=tie_tolerance,
        ),
        axis=0,
    )

    count_modal_index = np.argmax(selection_counts, axis=0)
    modal_index = count_modal_index.copy()
    count_tie = number_tied > 1
    resolvable = count_tie & np.isfinite(
        best_tie_score) & (score_tie_count == 1)
    modal_index[resolvable] = best_tie_index[resolvable]

    modal = np.take(
        np.asarray(ALL_PRIMARY_CODES, dtype=np.uint8), modal_index,
    ).astype(np.uint8)
    modal[count_tie & ~resolvable] = np.uint8(97)
    modal[available_count == 0] = NODATA_UINT8

    agreement = np.divide(
        100.0 * max_count,
        available_count,
        out=np.full(shape, np.nan, dtype=np.float64),
        where=available_count > 0,
    )
    return selection_counts, modal, agreement, available_count


def write_calibrated_archive(
    path: Path,
    source: Mapping[str, np.ndarray],
    winner: np.ndarray,
    member_available: np.ndarray,
    winning_margin: np.ndarray,
    calibrated_scores: Mapping[int, np.ndarray],
) -> None:
    """Write member-level calibrated scores without duplicating raw scores."""
    payload: dict[str, np.ndarray] = {
        'calibration_version': np.asarray(
            CALIBRATION_VERSION, dtype=np.int16,
        ),
        'calibration_method': np.asarray(CALIBRATION_METHOD),
        'calibration_zero_rule': np.asarray(CALIBRATION_ZERO_RULE),
        'primary_class': winner.astype(np.uint8),
        'available': member_available.astype(np.uint8),
        'connectivity_class': np.asarray(
            source['connectivity_class'], dtype=np.uint8,
        ),
        'winning_margin_calibrated_percentile_points': winning_margin.astype(
            np.float32,
        ),
    }
    for code, values in calibrated_scores.items():
        payload[f'calibrated_score_code_{code}'] = values.astype(
            np.float32, copy=False,
        )
    for key in (*PASSTHROUGH_ARRAY_KEYS, *MEMBER_DESCRIPTOR_KEYS):
        if key in source:
            payload[key] = np.asarray(source[key])

    temporary = path.with_name(path.stem + '.tmp.npz')
    np.savez_compressed(temporary, **payload)
    temporary.replace(path)


def process_group(
    group: Group,
    output_root: Path,
    calibration_manifest: Mapping,
    references: Mapping[int, np.ndarray],
    expected_members: int | None,
    tie_tolerance: float,
    chunk_cells: int,
    preserve_mode: str,
    write_member_maps: bool,
    write_archive: bool,
    overwrite_output: bool,
) -> None:
    """Recalculate all winner-dependent products for one SSP-by-radius group."""
    final_dir = output_root / group.name
    if final_dir.exists():
        if not overwrite_output:
            raise FileExistsError(
                f'{final_dir} exists. Use --overwrite-output to replace it.'
            )
        shutil.rmtree(final_dir)

    temporary_dir = output_root / f'.{group.name}.tmp-{os.getpid()}'
    if temporary_dir.exists():
        shutil.rmtree(temporary_dir)
    temporary_dir.mkdir(parents=True)

    try:
        preserved = preserve_existing_products(
            group.directory, temporary_dir, preserve_mode,
        )
        (temporary_dir / 'raw_archive_source.txt').write_text(
            str(group.archive.resolve()) + '\n', encoding='utf-8',
        )

        with np.load(group.archive, allow_pickle=False) as source:
            n_members, height, width = validate_archive_structure(
                source, group, expected_members,
            )
            shape = (n_members, height, width)
            spatial_shape = (height, width)
            available = np.asarray(source['available'], dtype=bool)
            connectivity = np.asarray(source['connectivity_class'])
            profile = template_profile(
                find_template_raster(group), spatial_shape,
            )
            inventory = load_inventory(group, source, n_members)

            best = np.full(shape, -np.inf, dtype=np.float32)
            second = np.full(shape, -np.inf, dtype=np.float32)
            best_code = np.zeros(shape, dtype=np.uint8)
            any_calibrated = np.zeros(shape, dtype=bool)
            raw_means: dict[int, np.ndarray] = {}
            calibrated_means: dict[int, np.ndarray] = {}
            calibrated_scores: dict[int, np.ndarray] = {}

            for code in POLICY_CODES:
                print(
                    f'{group.name}: calibrating policy code {code}',
                    flush=True,
                )
                raw = np.asarray(
                    source[f'common_score_code_{code}'], dtype=np.float32,
                )
                eligible = eligible_score_mask(
                    code, raw, available, connectivity,
                )
                # Validation also confirms that only negligible numerical
                # overshoot is clipped before ECDF mapping.
                checked_raw_scores(
                    raw, eligible, f'{group.name}, policy code {code}',
                )
                calibrated = calibrate_scores(
                    raw, eligible, references[code], chunk_cells,
                )
                calibrated_scores[code] = calibrated
                raw_masked = np.where(eligible, raw, np.nan).astype(np.float32)
                raw_means[code] = member_mean(raw_masked)
                calibrated_means[code] = member_mean(calibrated)
                any_calibrated |= np.isfinite(calibrated)
                update_top_two(
                    best, second, best_code, calibrated, code,
                )

            member_available = available & any_calibrated
            winner_count = np.zeros(shape, dtype=np.uint8)
            for calibrated in calibrated_scores.values():
                winner_count += (
                    np.isfinite(calibrated)
                    & np.isclose(
                        calibrated,
                        best,
                        rtol=0.0,
                        atol=tie_tolerance,
                    )
                ).astype(np.uint8)

            winner = np.full(shape, NODATA_UINT8, dtype=np.uint8)
            positive = member_available & (best > tie_tolerance)
            winner[member_available & ~positive] = np.uint8(0)
            unique = positive & (winner_count == 1)
            winner[unique] = best_code[unique]
            winner[positive & (winner_count > 1)] = np.uint8(97)

            winning_margin = np.full(shape, np.nan, dtype=np.float32)
            margin_valid = member_available & np.isfinite(
                best) & np.isfinite(second)
            winning_margin[margin_valid] = (
                best[margin_valid] - second[margin_valid]
            )

            (
                selection_counts,
                modal_class,
                modal_agreement,
                available_count,
            ) = compute_modal_products(
                winner,
                member_available,
                calibrated_means,
                tie_tolerance,
            )
            mean_margin = member_mean(winning_margin)

            if write_member_maps:
                member_dir = temporary_dir / 'member_priority_maps'
                member_dir.mkdir(parents=True, exist_ok=True)
                updated_inventory: list[dict[str, str]] = []
                for index, row in enumerate(inventory):
                    filename = member_map_name(row, index)
                    path = member_dir / filename
                    write_uint8(
                        path,
                        winner[index],
                        profile,
                        'Ensemble-member KM-GBF priority after fixed '
                        'positive-only pooled target calibration',
                    )
                    updated = dict(row)
                    updated['priority_map'] = str(
                        final_dir / 'member_priority_maps' / filename
                    )
                    updated_inventory.append(updated)
            else:
                updated_inventory = [dict(row) for row in inventory]
                for row in updated_inventory:
                    row['priority_map'] = 'not written (--skip-member-maps)'

            save_inventory(
                temporary_dir / 'member_inventory.csv', updated_inventory,
            )

            write_uint8(
                temporary_dir / 'kmgbf_modal_priority.tif',
                modal_class,
                profile,
                f'Modal calibrated KM-GBF policy class across {n_members} '
                'ensemble members',
            )
            write_float(
                temporary_dir / 'kmgbf_modal_priority_agreement_percent.tif',
                modal_agreement,
                profile,
                'Percentage of available ensemble members supporting the '
                'modal calibrated KM-GBF class',
            )
            write_float(
                temporary_dir / 'available_member_count.tif',
                available_count.astype(np.float64),
                profile,
                'Number of members with at least one eligible finite '
                'calibrated target score',
            )
            write_float(
                temporary_dir
                / 'kmgbf_mean_winning_margin_fixed_pooled_calibrated.tif',
                mean_margin,
                profile,
                'Mean difference between the highest and second-highest fixed '
                'positive-only pooled calibrated target scores, in percentile '
                'points',
            )

            for code in POLICY_CODES:
                write_float(
                    temporary_dir
                    / f'policy_mean_raw_common_score_code{code}.tif',
                    raw_means[code],
                    profile,
                    f'Mean raw common-range ecological score for policy code '
                    f'{code}; retained as magnitude evidence and not used '
                    'directly for winner selection',
                )
                write_float(
                    temporary_dir
                    / (
                        'policy_mean_fixed_pooled_calibrated_score_'
                        f'code{code}.tif'
                    ),
                    calibrated_means[code],
                    profile,
                    f'Mean fixed positive-only pooled target-specific '
                    f'calibrated score for policy code {code}',
                )

            for code_index, code in enumerate(ALL_PRIMARY_CODES):
                frequency = np.divide(
                    100.0 * selection_counts[code_index],
                    available_count,
                    out=np.full(spatial_shape, np.nan, dtype=np.float64),
                    where=available_count > 0,
                )
                write_float(
                    temporary_dir / f'policy_selection_percent_code{code}.tif',
                    frequency,
                    profile,
                    f'Percentage of available members selecting calibrated '
                    f'policy code {code}',
                )

            if write_archive:
                write_calibrated_archive(
                    temporary_dir
                    / 'kmgbf_member_scores_fixed_pooled_calibrated.npz',
                    source,
                    winner,
                    member_available,
                    winning_margin,
                    calibrated_scores,
                )

            group_metadata = {
                'created_utc': utc_now(),
                'ssp': group.ssp,
                'radius_pixels': group.radius,
                'radius_km': group.radius * 50,
                'member_count': n_members,
                'source_group': str(group.directory.resolve()),
                'raw_member_score_archive': str(group.archive.resolve()),
                'raw_scores_preserved': True,
                'raw_score_role': (
                    'Ecological magnitude evidence on the original common '
                    '0-100 range; not used directly to select the winner.'
                ),
                'calibration_manifest': str(
                    (output_root / 'calibration' / 'calibration_manifest.json')
                    .resolve()
                ),
                'calibration_method': calibration_manifest['method'],
                'calibration_version': calibration_manifest[
                    'calibration_version'
                ],
                'calibration_reference_population': calibration_manifest[
                    'reference_population'
                ],
                'calibration_zero_treatment': calibration_manifest['zero_rule'],
                'calibration_midrank_formula': calibration_manifest[
                    'midrank_formula'
                ],
                'calibration_target_audit_counts': {
                    code: {
                        'eligible_observation_count': record[
                            'eligible_observation_count'
                        ],
                        'eligible_zero_count': record['eligible_zero_count'],
                        'positive_reference_count': record[
                            'positive_reference_count'
                        ],
                    }
                    for code, record in calibration_manifest['targets'].items()
                },
                'same_calibration_for_all_ssps_and_radii': True,
                'winner_source': (
                    'maximum fixed positive-only pooled target-specific '
                    'calibrated score'
                ),
                'positive_tie_rule': (
                    f'Scores within absolute tolerance {tie_tolerance} of the '
                    'member maximum are tied; unresolved positive ties = 97.'
                ),
                'no_positive_action_code': 0,
                'target2_policy_code': 21,
                'unresolved_positive_tie_code': 97,
                'hotspot_percentiles': (
                    'Retained from the source group as within-group diagnostics '
                    'only; never used in calibrated winner selection.'
                ),
                'recipient_and_recovery_diagnostics': (
                    'Retained unchanged from source outputs. Source recovery '
                    'arrays are copied into the calibrated member archive when '
                    'present. They do not determine the winner.'
                ),
                'recalculated_from_calibrated_scores': [
                    'member priority maps',
                    'member winning margins',
                    'policy selection percentages',
                    'modal priority',
                    'modal priority agreement',
                    'mean calibrated target scores',
                ],
                'preserved_source_products': preserved,
                'calibrated_archive_written': write_archive,
                'member_maps_written': write_member_maps,
                'policy_labels': POLICY_LABELS,
            }
            atomic_write_json(
                temporary_dir / 'metadata.json', group_metadata,
            )

            # Free the largest arrays before moving the completed group.
            calibrated_scores.clear()

        temporary_dir.replace(final_dir)
        print(f'COMPLETED {group.name} -> {final_dir}', flush=True)

    except Exception:
        print(
            f'FAILED {group.name}; incomplete work retained at {temporary_dir}',
            file=sys.stderr,
            flush=True,
        )
        raise


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Command-line interface."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('--source-root', type=Path,
                        default=SOURCE_ROOT)
    parser.add_argument('--output-root', type=Path,
                        default=OUTPUT_ROOT)
    parser.add_argument(
        '--calibration-ssp', nargs='+', default=list(SSPS),
        help='SSPs pooled to build the fixed reference distributions.',
    )
    parser.add_argument(
        '--calibration-radius', nargs='+', type=int,
        default=list(RADII),
        help='Radii pooled to build the fixed reference distributions.',
    )
    parser.add_argument(
        '--apply-ssp', nargs='+', default=list(SSPS),
        help='SSPs to which the frozen calibration is applied.',
    )
    parser.add_argument(
        '--apply-radius', nargs='+', type=int, default=list(RADII),
        help='Radii to which the frozen calibration is applied.',
    )
    parser.add_argument(
        '--expected-members', type=int, default=180,
        help='Set to 0 to disable the member-count check.',
    )
    parser.add_argument(
        '--tie-tolerance', type=float, default=TIE_TOLERANCE,
    )
    parser.add_argument(
        '--chunk-cells', type=int, default=CHUNK_CELLS,
        help='Maximum flattened observations calibrated in one search chunk.',
    )
    parser.add_argument(
        '--preserve-mode', choices=('copy', 'hardlink', 'symlink', 'none'),
        default='copy',
        help='How unchanged source diagnostics are retained in the new root.',
    )
    parser.add_argument(
        '--reuse-calibration', action='store_true',
        help=(
            'Use the existing frozen positive-only calibration instead of '
            'rebuilding it.'
        ),
    )
    parser.add_argument(
        '--overwrite-calibration', action='store_true',
        help='Explicitly replace an existing calibration manifest and arrays.',
    )
    parser.add_argument(
        '--build-only', action='store_true',
        help='Build and validate the fixed calibration, but do not apply it.',
    )
    parser.add_argument(
        '--overwrite-output', action='store_true',
        help='Replace existing calibrated SSP-by-radius output groups.',
    )
    parser.add_argument(
        '--skip-member-maps', action='store_true',
        help='Diagnostic option: do not write the recalculated member GeoTIFFs.',
    )
    parser.add_argument(
        '--skip-calibrated-archive', action='store_true',
        help='Diagnostic option: do not write member-level calibrated scores.',
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    if args.tie_tolerance < 0:
        raise ValueError('--tie-tolerance must be non-negative')
    if args.chunk_cells < 1:
        raise ValueError('--chunk-cells must be positive')
    expected_members = args.expected_members or None

    source_root = args.source_root.resolve()
    output_root = args.output_root.resolve()
    if output_root == source_root:
        raise ValueError(
            '--output-root must differ from --source-root so raw products are '
            'never overwritten'
        )
    calibration_dir = output_root / 'calibration'
    output_root.mkdir(parents=True, exist_ok=True)

    calibration_groups = discover_groups(
        source_root,
        args.calibration_ssp,
        args.calibration_radius,
        expected_members,
    )

    if args.reuse_calibration:
        if args.overwrite_calibration:
            raise ValueError(
                '--reuse-calibration and --overwrite-calibration are mutually '
                'exclusive'
            )
        manifest, references = load_fixed_calibrations(calibration_dir)
        print(
            f'REUSING fixed positive-only calibration {calibration_dir}',
            flush=True,
        )
    else:
        manifest = build_fixed_calibrations(
            calibration_groups,
            calibration_dir,
            overwrite=args.overwrite_calibration,
        )
        manifest, references = load_fixed_calibrations(calibration_dir)

    if args.build_only:
        print('CALIBRATION BUILT AND VALIDATED; no groups applied', flush=True)
        return

    apply_groups = discover_groups(
        source_root,
        args.apply_ssp,
        args.apply_radius,
        expected_members,
    )
    for group in apply_groups:
        process_group(
            group=group,
            output_root=output_root,
            calibration_manifest=manifest,
            references=references,
            expected_members=expected_members,
            tie_tolerance=args.tie_tolerance,
            chunk_cells=args.chunk_cells,
            preserve_mode=args.preserve_mode,
            write_member_maps=not args.skip_member_maps,
            write_archive=not args.skip_calibrated_archive,
            overwrite_output=args.overwrite_output,
        )

    print(
        'ALL REQUESTED FIXED POSITIVE-ONLY POOLED TARGET CALIBRATIONS COMPLETED',
        flush=True,
    )


if __name__ == '__main__':
    main()
