#!/usr/bin/env python3
"""
Validate outputs from the fixed positive-only pooled target-specific
KM-GBF calibration workflow.

The validator checks:

1. Every requested SSP-by-radius output directory exists.
2. Required calibrated products are present.
3. The calibration manifest and positive-only sorted reference arrays are valid.
4. Group metadata agrees with the requested SSP, radius, calibration version,
   member count and calibration method.
5. All principal GeoTIFF products use the same grid.
6. Modal priority values use only documented policy codes.
7. Available-member counts are between zero and the expected ensemble size.
8. Percent rasters stay within 0-100, allowing a small numerical tolerance.
9. Raw mean scores and calibrated mean scores stay within 0-100.
10. Policy-selection percentages sum to approximately 100% where members are
    available.
11. Modal agreement equals the maximum policy-selection percentage.
12. The calibrated member archive contains the expected arrays and dimensions.
13. The member-level primary classes reproduce the modal priority and modal
    agreement products.
14. Eligible raw zeros and positive-only reference counts are internally
    consistent in the calibration manifest.
15. Optional member GeoTIFFs agree with the calibrated member archive.

This validator does not require source diagnostics that were merely preserved
from the earlier workflow because their presence depends on --preserve-mode.
Such diagnostics can be added to OPTIONAL_PRESERVED_PRODUCTS if required.
"""

from __future__ import annotations

import argparse
import csv
import json
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
    0: "No positive action",
    1: "T1 spatial planning",
    21: "T2 restoration",
    8: "T8 climate adaptation",
    10: "T10 production landscapes",
    31: "T31 maintain diffuse landscapes",
    32: "T32 safeguard concentrated connectivity",
    33: "T33 support PA/OECM under climate vulnerability",
    34: "T34 safeguard PA/OECM under land-use pressure",
    97: "Unresolved positive tie",
}

SSPS = ("ssp126", "ssp370", "ssp585")
RADII = (6, 14, 20)

ANALYSIS_ROOT = Path(
    "/capacity/occr_davin/mguzman/P1/Connectivity/approach2/"
    "kmgbf_memberwise_t_calibration_positive_180"
)

EXPECTED_MEMBER_COUNT = 180
EXPECTED_CALIBRATION_VERSION = 2

EXPECTED_CALIBRATION_METHOD = (
    "exact positive-only pooled target-specific empirical mid-rank ECDF"
)

NODATA_UINT8 = 255
PERCENT_TOLERANCE = 1.0e-3
SUM_TOLERANCE = 5.0e-2
AGREEMENT_TOLERANCE = 5.0e-2
SORT_CHECK_CHUNK = 1_000_000


# -----------------------------------------------------------------------------
# Required outputs
# -----------------------------------------------------------------------------

REQUIRED_GROUP_FILES = (
    "kmgbf_modal_priority.tif",
    "kmgbf_modal_priority_agreement_percent.tif",
    "available_member_count.tif",
    "kmgbf_mean_winning_margin_fixed_pooled_calibrated.tif",
    "kmgbf_member_scores_fixed_pooled_calibrated.npz",
    "member_inventory.csv",
    "metadata.json",
    "raw_archive_source.txt",
)

# These are expected only if they were present in the source group and the
# calibrated run used a preserve mode other than "none". They are therefore
# reported but not treated as unconditional failures.
OPTIONAL_PRESERVED_PRODUCTS = (
    "connectivity_impeded_agreement_percent.tif",
    "connectivity_diffuse_agreement_percent.tif",
    "connectivity_channelized_agreement_percent.tif",
    "connectivity_intensified_agreement_percent.tif",
    "connectivity_concentrated_agreement_percent.tif",
    "target2_mean_recovery_score.tif",
    "target2_mean_recovery_alignment_index.tif",
    "target2_mean_restoration_gap_index.tif",
)

REQUIRED_ARCHIVE_KEYS = {
    "calibration_version",
    "calibration_method",
    "calibration_zero_rule",
    "primary_class",
    "available",
    "connectivity_class",
    "winning_margin_calibrated_percentile_points",
    *{f"calibrated_score_code_{code}" for code in POLICY_CODES},
}


# -----------------------------------------------------------------------------
# Utility functions
# -----------------------------------------------------------------------------

def add_failure(failures: list[str], message: str) -> None:
    """Store a validation failure."""
    failures.append(message)


def load_json(path: Path) -> dict:
    """Read a JSON object."""
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path}: expected a JSON object")
    return value


def raster_signature(path: Path) -> tuple:
    """Return the grid definition used for raster comparisons."""
    with rasterio.open(path) as source:
        return (
            source.crs,
            source.transform,
            source.width,
            source.height,
        )


def read_masked(path: Path) -> np.ma.MaskedArray:
    """Read one raster band as a masked array."""
    with rasterio.open(path) as source:
        if source.count != 1:
            raise ValueError(
                f"{path}: expected one raster band, found {source.count}"
            )
        return source.read(1, masked=True)


def finite_values(path: Path) -> np.ndarray:
    """Return all unmasked finite raster values."""
    values = read_masked(path).compressed()
    return values[np.isfinite(values)]


def check_grid(
    path: Path,
    expected_signature: tuple,
    failures: list[str],
) -> None:
    """Check that a raster uses the expected grid."""
    try:
        signature = raster_signature(path)
    except Exception as error:
        add_failure(failures, f"Cannot read raster grid {path}: {error}")
        return

    if signature != expected_signature:
        add_failure(failures, f"Grid mismatch: {path}")


def check_range(
    path: Path,
    minimum: float,
    maximum: float,
    failures: list[str],
    tolerance: float = PERCENT_TOLERANCE,
) -> None:
    """Check the finite range of a continuous raster."""
    try:
        values = finite_values(path)
    except Exception as error:
        add_failure(failures, f"Cannot read {path}: {error}")
        return

    if values.size == 0:
        add_failure(failures, f"No finite raster values in {path}")
        return

    observed_minimum = float(np.min(values))
    observed_maximum = float(np.max(values))

    if observed_minimum < minimum - tolerance:
        add_failure(
            failures,
            f"{path}: minimum {observed_minimum:.8g} is below {minimum}",
        )

    if observed_maximum > maximum + tolerance:
        add_failure(
            failures,
            f"{path}: maximum {observed_maximum:.8g} exceeds {maximum}",
        )


def check_nonnegative(
    path: Path,
    failures: list[str],
    tolerance: float = PERCENT_TOLERANCE,
) -> None:
    """Check that all finite values are non-negative."""
    try:
        values = finite_values(path)
    except Exception as error:
        add_failure(failures, f"Cannot read {path}: {error}")
        return

    if values.size == 0:
        add_failure(failures, f"No finite raster values in {path}")
        return

    minimum = float(np.min(values))
    if minimum < -tolerance:
        add_failure(
            failures,
            f"{path}: minimum {minimum:.8g} is negative",
        )


def masked_to_nan(values: np.ma.MaskedArray) -> np.ndarray:
    """Convert a masked raster to float64 with NaN outside valid cells."""
    return np.asarray(values.filled(np.nan), dtype=np.float64)


def scalar_text(value: np.ndarray) -> str:
    """Convert a scalar NumPy archive value to text."""
    array = np.asarray(value)
    if array.size != 1:
        raise ValueError(f"Expected a scalar value, found shape {array.shape}")
    return str(array.reshape(()).item())


def scalar_int(value: np.ndarray) -> int:
    """Convert a scalar NumPy archive value to an integer."""
    array = np.asarray(value)
    if array.size != 1:
        raise ValueError(f"Expected a scalar value, found shape {array.shape}")
    return int(array.reshape(()).item())


def expected_group_files() -> list[str]:
    """Construct all unconditional group-level output filenames."""
    files = list(REQUIRED_GROUP_FILES)

    for code in POLICY_CODES:
        files.extend(
            [
                f"policy_mean_raw_common_score_code{code}.tif",
                (
                    "policy_mean_fixed_pooled_calibrated_score_"
                    f"code{code}.tif"
                ),
            ]
        )

    for code in ALL_PRIMARY_CODES:
        files.append(f"policy_selection_percent_code{code}.tif")

    return files


def compute_modal_from_archive(
    primary_class: np.ndarray,
    available: np.ndarray,
    calibrated_means: Mapping[int, np.ndarray],
    tie_tolerance: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Reproduce the calibrated workflow's modal class and agreement calculation.

    Ties in member-selection counts are resolved with the mean calibrated score
    for the tied policy classes. Code 97 is retained when that secondary
    comparison does not identify one unique class.
    """
    spatial_shape = primary_class.shape[1:]

    available_count = np.sum(
        available,
        axis=0,
        dtype=np.uint16,
    )

    selection_counts = np.stack(
        [
            np.sum(
                (primary_class == code) & available,
                axis=0,
                dtype=np.uint16,
            )
            for code in ALL_PRIMARY_CODES
        ]
    )

    max_count = np.max(selection_counts, axis=0)
    tied = selection_counts == max_count[None, :, :]
    number_tied = np.sum(tied, axis=0)

    mean_stack = np.stack(
        [
            (
                np.nan_to_num(
                    calibrated_means[code],
                    nan=-np.inf,
                )
                if code in calibrated_means
                else np.full(
                    spatial_shape,
                    -np.inf,
                    dtype=np.float64,
                )
            )
            for code in ALL_PRIMARY_CODES
        ]
    )

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

    modal_index = np.argmax(selection_counts, axis=0)
    count_tie = number_tied > 1

    resolvable = (
        count_tie
        & np.isfinite(best_tie_score)
        & (score_tie_count == 1)
    )
    modal_index[resolvable] = best_tie_index[resolvable]

    modal = np.take(
        np.asarray(ALL_PRIMARY_CODES, dtype=np.uint8),
        modal_index,
    ).astype(np.uint8)

    modal[count_tie & ~resolvable] = np.uint8(97)
    modal[available_count == 0] = np.uint8(NODATA_UINT8)

    agreement = np.divide(
        100.0 * max_count,
        available_count,
        out=np.full(spatial_shape, np.nan, dtype=np.float64),
        where=available_count > 0,
    )

    return modal, agreement, available_count


# -----------------------------------------------------------------------------
# Calibration validation
# -----------------------------------------------------------------------------

def validate_calibration(
    root: Path,
    failures: list[str],
) -> dict | None:
    """Validate the frozen calibration manifest and reference arrays."""
    calibration_directory = root / "calibration"
    manifest_path = calibration_directory / "calibration_manifest.json"

    print("CHECK CALIBRATION")

    if not manifest_path.is_file():
        add_failure(failures, f"Missing {manifest_path}")
        return None

    try:
        manifest = load_json(manifest_path)
    except Exception as error:
        add_failure(failures, f"Cannot read {manifest_path}: {error}")
        return None

    version = manifest.get("calibration_version")
    method = manifest.get("method")

    print(f"  version: {version}")
    print(f"  method: {method}")

    if version != EXPECTED_CALIBRATION_VERSION:
        add_failure(
            failures,
            f"{manifest_path}: calibration_version {version} != "
            f"{EXPECTED_CALIBRATION_VERSION}",
        )

    if method != EXPECTED_CALIBRATION_METHOD:
        add_failure(
            failures,
            f"{manifest_path}: unexpected calibration method {method!r}",
        )

    if manifest.get("frozen_across_ssps_and_radii") is not True:
        add_failure(
            failures,
            f"{manifest_path}: calibration is not marked as frozen across "
            "SSPs and radii",
        )

    targets = manifest.get("targets")
    if not isinstance(targets, dict):
        add_failure(
            failures,
            f"{manifest_path}: targets must be a JSON object",
        )
        return manifest

    for code in POLICY_CODES:
        record = targets.get(str(code))
        if not isinstance(record, dict):
            add_failure(
                failures,
                f"{manifest_path}: missing target record for code {code}",
            )
            continue

        reference_name = record.get("reference_file")
        if not reference_name:
            add_failure(
                failures,
                f"{manifest_path}: no reference_file for policy code {code}",
            )
            continue

        reference_path = calibration_directory / str(reference_name)
        if not reference_path.is_file():
            add_failure(failures, f"Missing {reference_path}")
            continue

        eligible_count = int(
            record.get("eligible_observation_count", -1)
        )
        zero_count = int(record.get("eligible_zero_count", -1))
        positive_count = int(
            record.get("positive_reference_count", -1)
        )

        if eligible_count <= 0:
            add_failure(
                failures,
                f"{manifest_path}: policy code {code} has invalid eligible "
                f"count {eligible_count}",
            )

        if zero_count < 0:
            add_failure(
                failures,
                f"{manifest_path}: policy code {code} has invalid zero count "
                f"{zero_count}",
            )

        if positive_count <= 0:
            add_failure(
                failures,
                f"{manifest_path}: policy code {code} has invalid positive "
                f"reference count {positive_count}",
            )

        if eligible_count != zero_count + positive_count:
            add_failure(
                failures,
                f"{manifest_path}: policy code {code} has inconsistent audit "
                f"counts: eligible={eligible_count}, zero={zero_count}, "
                f"positive={positive_count}",
            )

        if record.get("reference_excludes_zeros") is not True:
            add_failure(
                failures,
                f"{manifest_path}: policy code {code} reference is not marked "
                "as zero-excluding",
            )

        try:
            reference = np.load(
                reference_path,
                mmap_mode="r",
                allow_pickle=False,
            )
        except Exception as error:
            add_failure(
                failures,
                f"Cannot read {reference_path}: {error}",
            )
            continue

        if reference.ndim != 1:
            add_failure(
                failures,
                f"{reference_path}: expected one dimension, found "
                f"{reference.shape}",
            )
            continue

        if reference.size != positive_count:
            add_failure(
                failures,
                f"{reference_path}: contains {reference.size} values; "
                f"manifest reports {positive_count}",
            )

        previous = -np.inf
        finite_count = 0
        nonpositive_count = 0
        unsorted = False

        for start in range(0, reference.size, SORT_CHECK_CHUNK):
            stop = min(start + SORT_CHECK_CHUNK, reference.size)
            block = np.asarray(reference[start:stop])

            finite = np.isfinite(block)
            finite_count += int(np.count_nonzero(finite))
            nonpositive_count += int(np.count_nonzero(block <= 0.0))

            if block.size:
                if float(block[0]) < previous:
                    unsorted = True
                if np.any(block[1:] < block[:-1]):
                    unsorted = True
                previous = float(block[-1])

        if finite_count != reference.size:
            add_failure(
                failures,
                f"{reference_path}: contains "
                f"{reference.size - finite_count} non-finite values",
            )

        if nonpositive_count:
            add_failure(
                failures,
                f"{reference_path}: contains {nonpositive_count} values "
                "less than or equal to zero",
            )

        if unsorted:
            add_failure(
                failures,
                f"{reference_path}: values are not sorted ascending",
            )

        if reference.size:
            observed_minimum = float(reference[0])
            observed_maximum = float(reference[-1])

            recorded_minimum = float(
                record.get("minimum_positive", np.nan)
            )
            recorded_maximum = float(
                record.get("maximum_positive", np.nan)
            )

            if not np.isclose(
                observed_minimum,
                recorded_minimum,
                rtol=0.0,
                atol=1.0e-6,
            ):
                add_failure(
                    failures,
                    f"{reference_path}: minimum {observed_minimum} does not "
                    f"match manifest {recorded_minimum}",
                )

            if not np.isclose(
                observed_maximum,
                recorded_maximum,
                rtol=0.0,
                atol=1.0e-6,
            ):
                add_failure(
                    failures,
                    f"{reference_path}: maximum {observed_maximum} does not "
                    f"match manifest {recorded_maximum}",
                )

        print(
            f"  code {code}: eligible={eligible_count:,} | "
            f"zeros={zero_count:,} | positive={positive_count:,}"
        )

    return manifest


# -----------------------------------------------------------------------------
# Group validation
# -----------------------------------------------------------------------------

def validate_metadata(
    path: Path,
    ssp: str,
    radius: int,
    expected_members: int,
    failures: list[str],
) -> dict | None:
    """Validate group metadata."""
    try:
        metadata = load_json(path)
    except Exception as error:
        add_failure(failures, f"Cannot read {path}: {error}")
        return None

    member_count = int(metadata.get("member_count", -1))
    metadata_ssp = metadata.get("ssp")
    metadata_radius = int(metadata.get("radius_pixels", -1))
    metadata_version = metadata.get("calibration_version")
    metadata_method = metadata.get("calibration_method")

    print(
        f"  metadata members: {member_count} | "
        f"version: {metadata_version}"
    )

    if metadata_ssp != ssp:
        add_failure(
            failures,
            f"{path}: SSP {metadata_ssp!r} != {ssp!r}",
        )

    if metadata_radius != radius:
        add_failure(
            failures,
            f"{path}: radius_pixels {metadata_radius} != {radius}",
        )

    if member_count != expected_members:
        add_failure(
            failures,
            f"{path}: member_count {member_count} != {expected_members}",
        )

    if metadata_version != EXPECTED_CALIBRATION_VERSION:
        add_failure(
            failures,
            f"{path}: calibration_version {metadata_version} != "
            f"{EXPECTED_CALIBRATION_VERSION}",
        )

    if metadata_method != EXPECTED_CALIBRATION_METHOD:
        add_failure(
            failures,
            f"{path}: unexpected calibration method "
            f"{metadata_method!r}",
        )

    if metadata.get("same_calibration_for_all_ssps_and_radii") is not True:
        add_failure(
            failures,
            f"{path}: same_calibration_for_all_ssps_and_radii is not true",
        )

    if metadata.get("target2_policy_code") != 21:
        add_failure(
            failures,
            f"{path}: target2_policy_code is not 21",
        )

    if metadata.get("unresolved_positive_tie_code") != 97:
        add_failure(
            failures,
            f"{path}: unresolved_positive_tie_code is not 97",
        )

    if metadata.get("no_positive_action_code") != 0:
        add_failure(
            failures,
            f"{path}: no_positive_action_code is not 0",
        )

    return metadata


def validate_inventory(
    path: Path,
    expected_members: int,
    expect_member_maps: bool,
    failures: list[str],
) -> list[dict[str, str]]:
    """Validate the recalibrated member inventory."""
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
    except Exception as error:
        add_failure(failures, f"Cannot read {path}: {error}")
        return []

    if len(rows) != expected_members:
        add_failure(
            failures,
            f"{path}: expected {expected_members} rows, found {len(rows)}",
        )

    required_columns = {
        "member_index",
        "member_id",
        "ssp",
        "radius_pixels",
        "connectivity_taxon",
        "priority_map",
    }

    if rows:
        missing_columns = required_columns.difference(rows[0])
        if missing_columns:
            add_failure(
                failures,
                f"{path}: missing columns {sorted(missing_columns)}",
            )

    if expect_member_maps:
        for index, row in enumerate(rows):
            priority_map = row.get("priority_map", "")
            if not priority_map:
                add_failure(
                    failures,
                    f"{path}: row {index} has no priority_map",
                )

    return rows


def validate_archive(
    archive_path: Path,
    expected_members: int,
    expected_spatial_shape: tuple[int, int],
    tie_tolerance: float,
    failures: list[str],
) -> dict | None:
    """
    Validate the calibrated archive and derive independently reconstructed
    modal products for comparison with the written rasters.
    """
    try:
        source_context = np.load(archive_path, allow_pickle=False)
    except Exception as error:
        add_failure(failures, f"Cannot open {archive_path}: {error}")
        return None

    with source_context as source:
        missing = sorted(REQUIRED_ARCHIVE_KEYS.difference(source.files))
        if missing:
            add_failure(
                failures,
                f"{archive_path}: missing arrays {missing}",
            )
            return None

        try:
            version = scalar_int(source["calibration_version"])
            method = scalar_text(source["calibration_method"])
        except Exception as error:
            add_failure(
                failures,
                f"{archive_path}: invalid calibration metadata: {error}",
            )
            return None

        if version != EXPECTED_CALIBRATION_VERSION:
            add_failure(
                failures,
                f"{archive_path}: calibration_version {version} != "
                f"{EXPECTED_CALIBRATION_VERSION}",
            )

        if method != EXPECTED_CALIBRATION_METHOD:
            add_failure(
                failures,
                f"{archive_path}: unexpected calibration method {method!r}",
            )

        primary_class = np.asarray(source["primary_class"])
        available = np.asarray(source["available"], dtype=bool)
        connectivity = np.asarray(source["connectivity_class"])
        winning_margin = np.asarray(
            source["winning_margin_calibrated_percentile_points"],
            dtype=np.float32,
        )

        expected_shape = (
            expected_members,
            expected_spatial_shape[0],
            expected_spatial_shape[1],
        )

        for key, values in (
            ("primary_class", primary_class),
            ("available", available),
            ("connectivity_class", connectivity),
            (
                "winning_margin_calibrated_percentile_points",
                winning_margin,
            ),
        ):
            if values.shape != expected_shape:
                add_failure(
                    failures,
                    f"{archive_path}: {key} has shape {values.shape}; "
                    f"expected {expected_shape}",
                )

        if primary_class.shape != expected_shape:
            return None

        allowed = set(ALL_PRIMARY_CODES)
        observed_classes = set(
            np.unique(primary_class[available]).astype(int).tolist()
        )
        unexpected = observed_classes.difference(allowed)
        if unexpected:
            add_failure(
                failures,
                f"{archive_path}: unexpected primary classes "
                f"{sorted(unexpected)}",
            )

        unavailable_values = np.unique(primary_class[~available])
        invalid_unavailable = unavailable_values[
            unavailable_values != NODATA_UINT8
        ]
        if invalid_unavailable.size:
            add_failure(
                failures,
                f"{archive_path}: unavailable cells contain non-nodata "
                f"classes {invalid_unavailable.astype(int).tolist()}",
            )

        finite_margins = winning_margin[np.isfinite(winning_margin)]
        if finite_margins.size and np.min(finite_margins) < -PERCENT_TOLERANCE:
            add_failure(
                failures,
                f"{archive_path}: winning margin contains negative values",
            )

        calibrated_scores: dict[int, np.ndarray] = {}
        calibrated_means: dict[int, np.ndarray] = {}

        for code in POLICY_CODES:
            key = f"calibrated_score_code_{code}"
            values = np.asarray(source[key], dtype=np.float32)

            if values.shape != expected_shape:
                add_failure(
                    failures,
                    f"{archive_path}: {key} has shape {values.shape}; "
                    f"expected {expected_shape}",
                )
                continue

            finite = values[np.isfinite(values)]
            if finite.size:
                minimum = float(np.min(finite))
                maximum = float(np.max(finite))

                if minimum < -PERCENT_TOLERANCE:
                    add_failure(
                        failures,
                        f"{archive_path}: {key} minimum {minimum} is below 0",
                    )

                if maximum > 100.0 + PERCENT_TOLERANCE:
                    add_failure(
                        failures,
                        f"{archive_path}: {key} maximum {maximum} exceeds 100",
                    )

            calibrated_scores[code] = values

            finite_mask = np.isfinite(values)
            count = np.sum(finite_mask, axis=0)
            calibrated_means[code] = np.divide(
                np.nansum(values, axis=0, dtype=np.float64),
                count,
                out=np.full(
                    expected_spatial_shape,
                    np.nan,
                    dtype=np.float64,
                ),
                where=count > 0,
            )

        if len(calibrated_scores) != len(POLICY_CODES):
            return None

        any_score = np.zeros(expected_shape, dtype=bool)
        for values in calibrated_scores.values():
            any_score |= np.isfinite(values)

        if not np.array_equal(available, available & any_score):
            bad_count = int(np.count_nonzero(available & ~any_score))
            add_failure(
                failures,
                f"{archive_path}: {bad_count} available cells have no finite "
                "calibrated score",
            )

        modal, agreement, available_count = compute_modal_from_archive(
            primary_class,
            available,
            calibrated_means,
            tie_tolerance,
        )

        return {
            "primary_class": primary_class,
            "available": available,
            "calibrated_means": calibrated_means,
            "modal": modal,
            "agreement": agreement,
            "available_count": available_count,
        }


def validate_selection_percentages(
    folder: Path,
    available_count: np.ndarray,
    modal_agreement: np.ndarray,
    failures: list[str],
) -> None:
    """Validate selection-percentage sums and modal agreement."""
    selection_stack = []

    for code in ALL_PRIMARY_CODES:
        path = folder / f"policy_selection_percent_code{code}.tif"
        if not path.is_file():
            return

        try:
            values = masked_to_nan(read_masked(path))
        except Exception as error:
            add_failure(failures, f"Cannot read {path}: {error}")
            return

        selection_stack.append(values)

    stack = np.stack(selection_stack)
    valid = available_count > 0

    sums = np.nansum(stack, axis=0)
    bad_sum = valid & ~np.isclose(
        sums,
        100.0,
        rtol=0.0,
        atol=SUM_TOLERANCE,
    )

    if np.any(bad_sum):
        count = int(np.count_nonzero(bad_sum))
        maximum_error = float(np.max(np.abs(sums[bad_sum] - 100.0)))
        add_failure(
            failures,
            f"{folder}: policy-selection percentages do not sum to 100 in "
            f"{count} available cells; maximum error={maximum_error:.8g}",
        )

    written_maximum = np.nanmax(stack, axis=0)
    comparable = valid & np.isfinite(modal_agreement)

    mismatch = comparable & ~np.isclose(
        written_maximum,
        modal_agreement,
        rtol=0.0,
        atol=AGREEMENT_TOLERANCE,
    )

    if np.any(mismatch):
        count = int(np.count_nonzero(mismatch))
        maximum_error = float(
            np.max(
                np.abs(
                    written_maximum[mismatch]
                    - modal_agreement[mismatch]
                )
            )
        )
        add_failure(
            failures,
            f"{folder}: modal agreement differs from maximum policy-selection "
            f"percentage in {count} cells; maximum error="
            f"{maximum_error:.8g}",
        )


def validate_member_maps(
    folder: Path,
    inventory: Sequence[Mapping[str, str]],
    archive_primary_class: np.ndarray,
    reference_signature: tuple,
    failures: list[str],
) -> None:
    """Check optional member maps against the calibrated archive."""
    member_directory = folder / "member_priority_maps"

    if not member_directory.is_dir():
        add_failure(
            failures,
            f"Missing member-map directory {member_directory}",
        )
        return

    tif_paths = sorted(member_directory.glob("*.tif"))

    if len(tif_paths) != archive_primary_class.shape[0]:
        add_failure(
            failures,
            f"{member_directory}: expected "
            f"{archive_primary_class.shape[0]} GeoTIFFs, found "
            f"{len(tif_paths)}",
        )

    if len(inventory) != archive_primary_class.shape[0]:
        return

    for index, row in enumerate(inventory):
        inventory_path = row.get("priority_map", "")
        if not inventory_path:
            add_failure(
                failures,
                f"{folder / 'member_inventory.csv'}: row {index} has no "
                "priority_map",
            )
            continue

        path = member_directory / Path(inventory_path).name
        if not path.is_file():
            add_failure(failures, f"Missing {path}")
            continue

        check_grid(path, reference_signature, failures)

        try:
            written = np.asarray(read_masked(path).filled(NODATA_UINT8))
        except Exception as error:
            add_failure(failures, f"Cannot read {path}: {error}")
            continue

        expected = archive_primary_class[index]

        if not np.array_equal(written.astype(np.uint8), expected.astype(np.uint8)):
            difference_count = int(np.count_nonzero(written != expected))
            add_failure(
                failures,
                f"{path}: differs from calibrated archive member {index} in "
                f"{difference_count} cells",
            )


def validate_group(
    root: Path,
    ssp: str,
    radius: int,
    expected_members: int,
    tie_tolerance: float,
    check_member_maps: bool,
    failures: list[str],
) -> None:
    """Validate one calibrated SSP-by-radius output group."""
    folder = root / f"{ssp}_rad{radius}"
    print(f"\nCHECK {ssp} rad{radius}")

    if not folder.is_dir():
        add_failure(failures, f"Missing directory {folder}")
        return

    required_files = expected_group_files()
    missing = [
        folder / name
        for name in required_files
        if not (folder / name).is_file()
    ]

    for path in missing:
        add_failure(failures, f"Missing {path}")

    if missing:
        print(f"  missing required files: {len(missing)}")

    optional_present = [
        name
        for name in OPTIONAL_PRESERVED_PRODUCTS
        if (folder / name).is_file()
    ]
    print(
        f"  optional preserved diagnostics present: "
        f"{len(optional_present)}/{len(OPTIONAL_PRESERVED_PRODUCTS)}"
    )

    metadata_path = folder / "metadata.json"
    if metadata_path.is_file():
        validate_metadata(
            metadata_path,
            ssp,
            radius,
            expected_members,
            failures,
        )

    modal_path = folder / "kmgbf_modal_priority.tif"
    if not modal_path.is_file():
        return

    try:
        reference_signature = raster_signature(modal_path)
        with rasterio.open(modal_path) as source:
            spatial_shape = (source.height, source.width)
            modal_written = np.asarray(
                source.read(1, masked=True).filled(NODATA_UINT8),
                dtype=np.uint8,
            )
    except Exception as error:
        add_failure(failures, f"Cannot read {modal_path}: {error}")
        return

    for name in required_files:
        path = folder / name
        if path.suffix.lower() == ".tif" and path.is_file():
            check_grid(path, reference_signature, failures)

    for name in OPTIONAL_PRESERVED_PRODUCTS:
        path = folder / name

        if not path.is_file():
            continue

        # Ensure same spatial grid.
        check_grid(path, reference_signature, failures)

        # Connectivity agreement products are percentages and must therefore
        # remain within the closed interval 0-100 wherever values are valid.
        if (
            name.startswith("connectivity_")
            and name.endswith("_agreement_percent.tif")
        ):
            check_range(
                path,
                minimum=0.0,
                maximum=100.0,
                failures=failures,
            )

    allowed = set(ALL_PRIMARY_CODES)
    valid_modal = modal_written[modal_written != NODATA_UINT8]
    observed = set(np.unique(valid_modal).astype(int).tolist())
    unexpected = observed.difference(allowed)

    if unexpected:
        add_failure(
            failures,
            f"{modal_path}: unexpected policy codes {sorted(unexpected)}",
        )

    if valid_modal.size:
        unique, counts = np.unique(valid_modal, return_counts=True)
        print(
            "  modal:",
            " ".join(
                f"{int(code)}={100.0 * count / valid_modal.size:.2f}%"
                for code, count in zip(unique, counts)
            ),
        )
    else:
        add_failure(failures, f"{modal_path}: contains no valid policy cells")

    # -------------------------------------------------------------------------
    # Continuous raster ranges
    # -------------------------------------------------------------------------

    percent_paths = [
        folder / "kmgbf_modal_priority_agreement_percent.tif",
        *[
            folder / f"policy_selection_percent_code{code}.tif"
            for code in ALL_PRIMARY_CODES
        ],
        *[
            folder
            / (
                "policy_mean_fixed_pooled_calibrated_score_"
                f"code{code}.tif"
            )
            for code in POLICY_CODES
        ],
    ]

    for path in percent_paths:
        if path.is_file():
            check_range(path, 0.0, 100.0, failures)

    for code in POLICY_CODES:
        path = folder / f"policy_mean_raw_common_score_code{code}.tif"
        if path.is_file():
            check_range(path, 0.0, 100.0, failures)

    margin_path = (
        folder
        / "kmgbf_mean_winning_margin_fixed_pooled_calibrated.tif"
    )
    if margin_path.is_file():
        check_range(margin_path, 0.0, 100.0, failures)

    count_path = folder / "available_member_count.tif"
    available_count_written = None

    if count_path.is_file():
        try:
            count_raster = masked_to_nan(read_masked(count_path))
            finite = count_raster[np.isfinite(count_raster)]

            if finite.size == 0:
                add_failure(
                    failures,
                    f"{count_path}: contains no finite values",
                )
            else:
                minimum = float(np.min(finite))
                maximum = float(np.max(finite))

                print(
                    f"  available members per cell: "
                    f"min={minimum:.0f} max={maximum:.0f}"
                )

                if minimum < -PERCENT_TOLERANCE:
                    add_failure(
                        failures,
                        f"{count_path}: minimum {minimum} is below zero",
                    )

                if maximum > expected_members + PERCENT_TOLERANCE:
                    add_failure(
                        failures,
                        f"{count_path}: maximum {maximum} exceeds "
                        f"{expected_members}",
                    )

                if not np.allclose(
                    finite,
                    np.rint(finite),
                    rtol=0.0,
                    atol=1.0e-6,
                ):
                    add_failure(
                        failures,
                        f"{count_path}: contains non-integer count values",
                    )

            available_count_written = count_raster
        except Exception as error:
            add_failure(failures, f"Cannot read {count_path}: {error}")

    # -------------------------------------------------------------------------
    # Inventory and calibrated archive
    # -------------------------------------------------------------------------

    inventory_path = folder / "member_inventory.csv"
    inventory: list[dict[str, str]] = []

    if inventory_path.is_file():
        inventory = validate_inventory(
            inventory_path,
            expected_members,
            check_member_maps,
            failures,
        )

    archive_path = (
        folder / "kmgbf_member_scores_fixed_pooled_calibrated.npz"
    )
    reconstructed = None

    if archive_path.is_file():
        reconstructed = validate_archive(
            archive_path,
            expected_members,
            spatial_shape,
            tie_tolerance,
            failures,
        )

    if reconstructed is None:
        return

    # -------------------------------------------------------------------------
    # Compare reconstructed products against written rasters
    # -------------------------------------------------------------------------

    modal_expected = reconstructed["modal"]
    modal_difference = modal_written != modal_expected

    if np.any(modal_difference):
        count = int(np.count_nonzero(modal_difference))
        add_failure(
            failures,
            f"{modal_path}: differs from the calibrated archive in "
            f"{count} cells",
        )

    agreement_path = (
        folder / "kmgbf_modal_priority_agreement_percent.tif"
    )
    modal_agreement_written = None

    if agreement_path.is_file():
        try:
            modal_agreement_written = masked_to_nan(
                read_masked(agreement_path)
            )
            modal_agreement_expected = reconstructed["agreement"]

            comparable = (
                np.isfinite(modal_agreement_written)
                & np.isfinite(modal_agreement_expected)
            )

            mismatch = comparable & ~np.isclose(
                modal_agreement_written,
                modal_agreement_expected,
                rtol=0.0,
                atol=AGREEMENT_TOLERANCE,
            )

            missing_written = (
                np.isfinite(modal_agreement_expected)
                & ~np.isfinite(modal_agreement_written)
            )

            unexpected_written = (
                ~np.isfinite(modal_agreement_expected)
                & np.isfinite(modal_agreement_written)
            )

            if np.any(mismatch | missing_written | unexpected_written):
                count = int(
                    np.count_nonzero(
                        mismatch | missing_written | unexpected_written
                    )
                )
                add_failure(
                    failures,
                    f"{agreement_path}: differs from the calibrated archive "
                    f"in {count} cells",
                )
        except Exception as error:
            add_failure(
                failures,
                f"Cannot compare {agreement_path}: {error}",
            )

    available_count_expected = reconstructed["available_count"]

    if available_count_written is not None:
        comparable = np.isfinite(available_count_written)

        mismatch = comparable & ~np.isclose(
            available_count_written,
            available_count_expected,
            rtol=0.0,
            atol=1.0e-6,
        )

        if np.any(mismatch):
            count = int(np.count_nonzero(mismatch))
            add_failure(
                failures,
                f"{count_path}: differs from the calibrated archive in "
                f"{count} cells",
            )

    if modal_agreement_written is not None:
        validate_selection_percentages(
            folder,
            available_count_expected,
            modal_agreement_written,
            failures,
        )

    # Compare written mean calibrated scores with archive-derived means.
    for code in POLICY_CODES:
        path = (
            folder
            / (
                "policy_mean_fixed_pooled_calibrated_score_"
                f"code{code}.tif"
            )
        )
        if not path.is_file():
            continue

        try:
            written = masked_to_nan(read_masked(path))
        except Exception as error:
            add_failure(failures, f"Cannot read {path}: {error}")
            continue

        expected = reconstructed["calibrated_means"][code]

        both_finite = np.isfinite(written) & np.isfinite(expected)
        mismatch = both_finite & ~np.isclose(
            written,
            expected,
            rtol=0.0,
            atol=1.0e-4,
        )

        missing_written = np.isfinite(expected) & ~np.isfinite(written)
        unexpected_written = ~np.isfinite(expected) & np.isfinite(written)

        difference = mismatch | missing_written | unexpected_written

        if np.any(difference):
            count = int(np.count_nonzero(difference))
            add_failure(
                failures,
                f"{path}: differs from the calibrated archive mean in "
                f"{count} cells",
            )

    if check_member_maps:
        validate_member_maps(
            folder,
            inventory,
            reconstructed["primary_class"],
            reference_signature,
            failures,
        )


# -----------------------------------------------------------------------------
# Command-line interface
# -----------------------------------------------------------------------------

def parse_args(
    argv: Sequence[str] | None = None,
) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--analysis-root",
        type=Path,
        default=ANALYSIS_ROOT,
        help="Root containing calibration/ and SSP-by-radius output groups.",
    )

    parser.add_argument(
        "--ssp",
        nargs="+",
        choices=SSPS,
        default=list(SSPS),
        help="SSP output groups to validate.",
    )

    parser.add_argument(
        "--radius",
        nargs="+",
        type=int,
        choices=RADII,
        default=list(RADII),
        help="Connectivity radii to validate.",
    )

    parser.add_argument(
        "--expected-members",
        type=int,
        default=EXPECTED_MEMBER_COUNT,
        help="Expected ensemble-member count.",
    )

    parser.add_argument(
        "--tie-tolerance",
        type=float,
        default=1.0e-6,
        help="Tolerance used to reconstruct calibrated winner ties.",
    )

    parser.add_argument(
        "--check-member-maps",
        action="store_true",
        help=(
            "Validate every recalibrated member GeoTIFF against the calibrated "
            "member archive. This is slower but provides the strictest check."
        ),
    )

    parser.add_argument(
        "--skip-calibration-check",
        action="store_true",
        help=(
            "Skip validation of the pooled positive-only calibration manifest "
            "and sorted reference arrays."
        ),
    )

    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    """Run all requested validation checks."""
    args = parse_args(argv)

    if args.expected_members < 1:
        raise ValueError("--expected-members must be at least one")

    if args.tie_tolerance < 0:
        raise ValueError("--tie-tolerance must be non-negative")

    root = args.analysis_root.resolve()
    if not root.is_dir():
        raise FileNotFoundError(root)

    failures: list[str] = []

    print("=" * 78)
    print("KM-GBF FIXED POSITIVE-ONLY POOLED CALIBRATION VALIDATION")
    print("=" * 78)
    print(f"Root: {root}")
    print(f"Expected members: {args.expected_members}")
    print(f"SSPs: {', '.join(args.ssp)}")
    print(f"Radii: {', '.join(str(value) for value in args.radius)}")
    print(f"Check member maps: {args.check_member_maps}")

    if not args.skip_calibration_check:
        validate_calibration(root, failures)

    for ssp in args.ssp:
        for radius in args.radius:
            validate_group(
                root=root,
                ssp=ssp,
                radius=radius,
                expected_members=args.expected_members,
                tie_tolerance=args.tie_tolerance,
                check_member_maps=args.check_member_maps,
                failures=failures,
            )

    print("\n" + "=" * 78)

    if failures:
        print(f"FAILED: {len(failures)} validation problem(s)")
        print("=" * 78)

        for index, failure in enumerate(failures, start=1):
            print(f"{index:03d}. {failure}")

        raise SystemExit(1)

    print("ALL CHECKS PASSED")
    print("=" * 78)


if __name__ == "__main__":
    main()
