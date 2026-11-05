#!/usr/bin/env python3
"""Build web TIFF catalog from finalised SSP x radius outputs."""
import json
import shutil
import subprocess
from config import (
    ANALYSIS_ROOT, DOCS, WEB_DATA_DIR, SSPS, RADII,
    SSP_LABELS, RADIUS_LABELS, ACTION_CODES, ACTION_LABELS, ACTION_COLORS,
    EXPECTED_MEMBER_COUNT,
)

WEB_DATA_DIR.mkdir(parents=True, exist_ok=True)


def export(source, destination):
    if not source.exists():
        print("Missing:", source)
        return False
    if shutil.which("gdal_translate"):
        subprocess.run([
            "gdal_translate", "-of", "COG", "-co", "COMPRESS=DEFLATE",
            "-co", "OVERVIEWS=AUTO", str(source), str(destination)
        ], check=True)
    else:
        shutil.copy2(source, destination)
    return True


def add(layers, source, filename, title, kind, group, ssp, radius, **extra):
    if export(source, WEB_DATA_DIR / filename):
        layers.append({
            "file": f"data/{filename}", "title": title, "kind": kind,
            "group": group, "ssp": ssp, "radius": radius, **extra
        })


def read_ensemble_info():
    """Read member count and connectivity taxa from every group's metadata."""
    counts = {}
    taxa = None
    design = None
    for ssp in SSPS:
        for radius in RADII:
            path = ANALYSIS_ROOT / f"{ssp}_rad{radius}" / "metadata.json"
            if not path.is_file():
                raise FileNotFoundError(path)
            meta = json.loads(path.read_text(encoding="utf-8"))
            counts[f"{ssp}_rad{radius}"] = int(meta["member_count"])
            taxa = meta.get("connectivity_taxa", taxa)
            design = meta.get("ensemble_design", design)

    sizes = set(counts.values())
    if sizes != {EXPECTED_MEMBER_COUNT}:
        raise ValueError(
            f"Expected {EXPECTED_MEMBER_COUNT} members in every group; "
            f"found {counts}"
        )
    return EXPECTED_MEMBER_COUNT, taxa, design


def main():
    ensemble_size, connectivity_taxa, ensemble_design = read_ensemble_info()
    layers = []

    for ssp in SSPS:
        for radius in RADII:
            root = ANALYSIS_ROOT / f"{ssp}_rad{radius}"
            add(layers, root / "kmgbf_modal_priority.tif",
                f"priority_{ssp}_rad{radius}.tif",
                "Modal strategic action", "categorical", "priority", ssp, radius)
            add(layers, root / "kmgbf_modal_priority_agreement_percent.tif",
                f"modal_agreement_{ssp}_rad{radius}.tif",
                "Modal strategic-action ensemble agreement (%)", "continuous",
                "priority", ssp, radius,
                min=0, max=100, palette="agreement", metric="modal_agreement")

            for code in ACTION_CODES:
                label = ACTION_LABELS[code]
                add(layers, root / f"policy_selection_percent_code{code}.tif",
                    f"selection_code{code}_{ssp}_rad{radius}.tif",
                    f"{label} — ensemble agreement (%)", "continuous",
                    "targets", ssp, radius,
                    action=code, min=0, max=100, palette=f"action_{code}",
                    metric="ensemble_agreement")
                add(layers, root / f"policy_mean_normalized_score_code{code}.tif",
                    f"score_code{code}_{ssp}_rad{radius}.tif",
                    f"{label} — mean normalised priority score", "continuous",
                    "target_scores", ssp, radius,
                    action=code, min=0, max=100, palette="plasma",
                    metric="mean_priority_score")

            for name in ["impeded", "diffuse", "channelized", "intensified",
                         "concentrated"]:
                add(layers, root / f"connectivity_{name}_agreement_percent.tif",
                    f"connectivity_{name}_{ssp}_rad{radius}.tif",
                    f"{name.capitalize()} connectivity — ensemble agreement (%)",
                    "continuous", "connectivity", ssp, radius,
                    min=0, max=100, palette="agreement",
                    metric="connectivity_agreement")

            add(layers, root / "target2_mean_recovery_score.tif",
                f"t2_recovery_{ssp}_rad{radius}.tif",
                "Projected secondary-land recovery score (0–1)", "continuous",
                "target2", ssp, radius,
                min=0, max=1, palette="ylgn", metric="recovery_score")
            add(layers, root / "target2_mean_recovery_alignment_index.tif",
                f"t2_alignment_{ssp}_rad{radius}.tif",
                "Restoration priority–recovery alignment", "continuous",
                "target2", ssp, radius,
                min=0, max=100, palette="ylgn", metric="recovery_alignment")
            add(layers, root / "target2_mean_restoration_gap_index.tif",
                f"t2_gap_{ssp}_rad{radius}.tif",
                "Restoration priority gap", "continuous", "target2", ssp, radius,
                min=0, max=100, palette="magma", metric="restoration_gap")

    catalog = {
        "classes": {str(k): v for k, v in ACTION_LABELS.items()},
        "colors": {str(k): v for k, v in ACTION_COLORS.items()},
        "sspLabels": SSP_LABELS,
        "radiusLabels": {str(k): v for k, v in RADIUS_LABELS.items()},
        "projection": "EPSG:8857",
        "ensembleSize": ensemble_size,
        "ensembleDesign": ensemble_design,
        "connectivityTaxa": connectivity_taxa,
        "layers": layers,
    }
    (DOCS / "catalog.json").write_text(
        json.dumps(catalog, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Wrote {len(layers)} layers | ensemble size {ensemble_size}")


if __name__ == "__main__":
    main()
