#!/usr/bin/env python3
"""
Build the web GeoTIFF catalogue from the final KM-GBF outputs produced from
positive-only target-specific calibration (across SSP-RCPs and radii).

The catalogue includes:
- modal strategic action;
- modal ensemble agreement;
- action-specific ensemble selection percentages;
- mean fixed-pooled calibrated scores;
- mean raw ecological scores;
- connectivity agreement diagnostics;
- Target 2 recovery diagnostics;
- available-member count;
- mean calibrated winning margin.

"""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

from config import (
    ANALYSIS_ROOT,
    DOCS,
    WEB_DATA_DIR,
    SSPS,
    RADII,
    SSP_LABELS,
    RADIUS_LABELS,
    ACTION_CODES,
    TIE_CODE,
    ACTION_LABELS,
    ACTION_COLORS,
    EXPECTED_MEMBER_COUNT,
)


EXPECTED_CALIBRATION_VERSION = 2

EXPECTED_CALIBRATION_METHOD = (
    "exact positive-only pooled target-specific empirical mid-rank ECDF"
)


def export(source: Path, destination: Path) -> bool:
    """
    Export a GeoTIFF to the web-data directory.

    If GDAL is available, the raster is converted to a Cloud Optimised
    GeoTIFF. Otherwise, the original GeoTIFF is copied unchanged.
    """
    if not source.is_file():
        print("Missing optional layer:", source)
        return False

    destination.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    if destination.exists():
        destination.unlink()

    if shutil.which("gdal_translate"):
        subprocess.run(
            [
                "gdal_translate",
                "-of",
                "COG",
                "-co",
                "COMPRESS=DEFLATE",
                "-co",
                "OVERVIEWS=AUTO",
                str(source),
                str(destination),
            ],
            check=True,
        )
    else:
        shutil.copy2(
            source,
            destination,
        )

    return True


def add(
    layers,
    source,
    filename,
    title,
    kind,
    group,
    ssp,
    radius,
    **extra,
):
    """Export one raster and append its information to the web catalogue."""
    destination = WEB_DATA_DIR / filename

    if not export(source, destination):
        return

    layer = {
        "file": f"data/{filename}",
        "title": title,
        "kind": kind,
        "group": group,
        "ssp": ssp,
        "radius": radius,
    }

    layer.update(extra)
    layers.append(layer)


def read_ensemble_info():
    """
    Validate metadata for every SSP x radius group.

    The web catalogue should only be generated when every requested group has
    the expected 180 ensemble members and uses the same fixed calibration.
    """
    counts = {}
    versions = {}
    methods = {}

    for ssp in SSPS:
        for radius in RADII:
            group_name = f"{ssp}_rad{radius}"
            path = (
                ANALYSIS_ROOT
                / group_name
                / "metadata.json"
            )

            if not path.is_file():
                raise FileNotFoundError(path)

            metadata = json.loads(
                path.read_text(
                    encoding="utf-8",
                )
            )

            if "member_count" not in metadata:
                raise KeyError(
                    f"{path}: member_count is missing"
                )

            if "calibration_version" not in metadata:
                raise KeyError(
                    f"{path}: calibration_version is missing"
                )

            if "calibration_method" not in metadata:
                raise KeyError(
                    f"{path}: calibration_method is missing"
                )

            member_count = int(
                metadata["member_count"]
            )

            calibration_version = int(
                metadata["calibration_version"]
            )

            calibration_method = str(
                metadata["calibration_method"]
            )

            counts[group_name] = member_count
            versions[group_name] = calibration_version
            methods[group_name] = calibration_method

            metadata_ssp = metadata.get("ssp")

            if metadata_ssp != ssp:
                raise ValueError(
                    f"{path}: metadata SSP "
                    f"{metadata_ssp!r} does not match "
                    f"{ssp!r}"
                )

            metadata_radius = int(
                metadata.get(
                    "radius_pixels",
                    -1,
                )
            )

            if metadata_radius != radius:
                raise ValueError(
                    f"{path}: metadata radius "
                    f"{metadata_radius} does not match "
                    f"{radius}"
                )

            same_calibration = metadata.get(
                "same_calibration_for_all_ssps_and_radii"
            )

            if same_calibration is not True:
                raise ValueError(
                    f"{path}: calibration is not marked "
                    "as common across SSPs and radii"
                )

    for group_name, count in counts.items():
        if count != EXPECTED_MEMBER_COUNT:
            raise ValueError(
                f"{group_name}: expected "
                f"{EXPECTED_MEMBER_COUNT} members, "
                f"found {count}"
            )

    for group_name, version in versions.items():
        if version != EXPECTED_CALIBRATION_VERSION:
            raise ValueError(
                f"{group_name}: expected calibration "
                f"version {EXPECTED_CALIBRATION_VERSION}, "
                f"found {version}"
            )

    for group_name, method in methods.items():
        if method != EXPECTED_CALIBRATION_METHOD:
            raise ValueError(
                f"{group_name}: unexpected calibration "
                f"method {method!r}"
            )

    return {
        "ensembleSize": EXPECTED_MEMBER_COUNT,
        "calibrationVersion": EXPECTED_CALIBRATION_VERSION,
        "calibrationMethod": EXPECTED_CALIBRATION_METHOD,
        "groupMemberCounts": counts,
    }


def add_priority_layers(
    layers,
    root,
    ssp,
    radius,
):
    """Add overall priority and agreement products."""

    add(
        layers,
        root / "kmgbf_modal_priority.tif",
        f"priority_{ssp}_rad{radius}.tif",
        "Modal strategic action",
        "categorical",
        "priority",
        ssp,
        radius,
        metric="modal_priority",
    )

    add(
        layers,
        root
        / "kmgbf_modal_priority_agreement_percent.tif",
        f"modal_agreement_{ssp}_rad{radius}.tif",
        "Modal strategic-action ensemble agreement (%)",
        "continuous",
        "priority",
        ssp,
        radius,
        min=0,
        max=100,
        palette="agreement",
        metric="modal_agreement",
        units="percent",
    )

    add(
        layers,
        root
        / "kmgbf_mean_winning_margin_fixed_pooled_calibrated.tif",
        f"winning_margin_{ssp}_rad{radius}.tif",
        "Mean calibrated winning margin",
        "continuous",
        "priority_diagnostics",
        ssp,
        radius,
        min=0,
        max=100,
        palette="plasma",
        metric="winning_margin",
        units="percentile points",
    )


def add_action_layers(
    layers,
    root,
    ssp,
    radius,
):
    """
    Add action-specific ensemble, calibrated-score and raw-score products.
    """

    for code in ACTION_CODES:
        label = ACTION_LABELS[code]

        selection_source = (
            root
            / f"policy_selection_percent_code{code}.tif"
        )

        selection_filename = (
            f"selection_code{code}_"
            f"{ssp}_rad{radius}.tif"
        )

        add(
            layers,
            selection_source,
            selection_filename,
            f"{label} - ensemble selection (%)",
            "continuous",
            "targets",
            ssp,
            radius,
            action=code,
            min=0,
            max=100,
            palette=f"action_{code}",
            metric="ensemble_agreement",
            units="percent",
            interpretation=(
                "Percentage of available ensemble members "
                "selecting this action."
            ),
        )

        calibrated_source = (
            root
            / (
                "policy_mean_fixed_pooled_calibrated_"
                f"score_code{code}.tif"
            )
        )

        calibrated_filename = (
            f"calibrated_score_code{code}_"
            f"{ssp}_rad{radius}.tif"
        )

        add(
            layers,
            calibrated_source,
            calibrated_filename,
            f"{label} - mean calibrated score",
            "continuous",
            "calibrated_target_scores",
            ssp,
            radius,
            action=code,
            min=0,
            max=100,
            palette="plasma",
            metric="mean_calibrated_score",
            units="empirical percentile",
            interpretation=(
                "Mean member-level score after fixed "
                "positive-only pooled target-specific "
                "empirical calibration."
            ),
        )


def add_connectivity_layers(
    layers,
    root,
    ssp,
    radius,
):
    """
    Add preserved connectivity agreement products.

    """

    connectivity_layers = [
        (
            "impeded",
            "Impeded connectivity",
        ),
        (
            "diffuse",
            "Diffuse connectivity",
        ),
        (
            "channelized",
            "Channelised connectivity",
        ),
        (
            "intensified",
            "Intensified connectivity",
        ),
        (
            "concentrated",
            "Concentrated connectivity",
        ),
    ]

    for name, label in connectivity_layers:
        source = (
            root
            / f"connectivity_{name}_agreement_percent.tif"
        )

        filename = (
            f"connectivity_{name}_"
            f"{ssp}_rad{radius}.tif"
        )

        add(
            layers,
            source,
            filename,
            f"{label} - ensemble agreement (%)",
            "continuous",
            "connectivity",
            ssp,
            radius,
            min=0,
            max=100,
            palette="agreement",
            metric="connectivity_agreement",
            units="percent",
        )


def add_target2_layers(
    layers,
    root,
    ssp,
    radius,
):
    """
    Add retained Target 2 recovery diagnostics.

    These diagnostics are retained for interpretation but do not determine
    the calibrated winning action.
    """

    add(
        layers,
        root / "target2_mean_recovery_score.tif",
        f"t2_recovery_{ssp}_rad{radius}.tif",
        "Projected secondary-land recovery score (0-1)",
        "continuous",
        "target2",
        ssp,
        radius,
        min=0,
        max=1,
        palette="ylgn",
        metric="recovery_score",
        units="index",
    )

    add(
        layers,
        root
        / "target2_mean_recovery_alignment_index.tif",
        f"t2_alignment_{ssp}_rad{radius}.tif",
        "Restoration priority-recovery alignment",
        "continuous",
        "target2",
        ssp,
        radius,
        min=0,
        max=100,
        palette="ylgn",
        metric="recovery_alignment",
        units="index",
    )

    add(
        layers,
        root
        / "target2_mean_restoration_gap_index.tif",
        f"t2_gap_{ssp}_rad{radius}.tif",
        "Restoration priority gap",
        "continuous",
        "target2",
        ssp,
        radius,
        min=0,
        max=100,
        palette="magma",
        metric="restoration_gap",
        units="index",
    )


def main():
    """Generate all web GeoTIFFs and write catalog.json."""

    if not ANALYSIS_ROOT.is_dir():
        raise FileNotFoundError(
            ANALYSIS_ROOT
        )

    DOCS.mkdir(
        parents=True,
        exist_ok=True,
    )

    # Rebuild the web TIFF directory so obsolete products from the earlier
    # normalised-score workflow cannot remain in the online catalogue.
    if WEB_DATA_DIR.exists():
        shutil.rmtree(
            WEB_DATA_DIR
        )

    WEB_DATA_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    ensemble_info = read_ensemble_info()

    layers = []

    for ssp in SSPS:
        for radius in RADII:
            root = (
                ANALYSIS_ROOT
                / f"{ssp}_rad{radius}"
            )

            if not root.is_dir():
                raise FileNotFoundError(
                    root
                )

            print(
                f"BUILD {ssp} rad{radius}"
            )

            add_priority_layers(
                layers,
                root,
                ssp,
                radius,
            )

            add_action_layers(
                layers,
                root,
                ssp,
                radius,
            )

            add_connectivity_layers(
                layers,
                root,
                ssp,
                radius,
            )

            add_target2_layers(
                layers,
                root,
                ssp,
                radius,
            )

    classes = {}

    for code, label in ACTION_LABELS.items():
        classes[str(code)] = label

    colors = {}

    for code, color in ACTION_COLORS.items():
        colors[str(code)] = color

    radius_labels = {}

    for radius, label in RADIUS_LABELS.items():
        radius_labels[str(radius)] = label

    catalog = {
        "classes": classes,
        "colors": colors,
        "sspLabels": SSP_LABELS,
        "radiusLabels": radius_labels,
        "projection": "EPSG:8857",
        "ensembleSize": ensemble_info["ensembleSize"],
        "calibrationVersion": ensemble_info[
            "calibrationVersion"
        ],
        "calibrationMethod": ensemble_info[
            "calibrationMethod"
        ],
        "groupMemberCounts": ensemble_info[
            "groupMemberCounts"
        ],
        "scoreInterpretation": {
            "raw": (
                "Raw common-range ecological magnitude on a "
                "0-100 scale. Raw scores are retained as "
                "ecological evidence and are not directly used "
                "for cross-target winner selection."
            ),
            "calibrated": (
                "Fixed positive-only pooled target-specific "
                "empirical percentile. These calibrated scores "
                "are used for member-level winner selection."
            ),
            "selection": (
                "Percentage of available ensemble members for "
                "which the action is the selected member-level "
                "outcome."
            ),
        },
        "layers": layers,
    }

    catalog_path = (
        DOCS
        / "catalog.json"
    )

    catalog_path.write_text(
        json.dumps(
            catalog,
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )

    print()
    print(
        f"Wrote {len(layers)} layers"
    )
    print(
        f"Catalogue: {catalog_path}"
    )
    print(
        f"Ensemble size: "
        f"{ensemble_info['ensembleSize']}"
    )
    print(
        "Calibration: "
        f"{ensemble_info['calibrationMethod']}"
    )


if __name__ == "__main__":
    main()
