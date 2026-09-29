#!/usr/bin/env python3
"""Convert catalogued EPSG:8857 GeoTIFFs to transparent RGBA PNG overlays."""

from pathlib import Path
import json

import numpy as np
import rasterio
from PIL import Image
from matplotlib import colormaps

DOCS = Path("docs")
CATALOG_PATH = DOCS / "catalog.json"
PNG_DIR = DOCS / "data_png"
PNG_CATALOG_PATH = DOCS / "catalog_png.json"


def hex_rgba(value):
    value = value.lstrip("#")
    if len(value) == 6:
        value += "FF"
    return tuple(int(value[i:i + 2], 16) for i in (0, 2, 4, 6))


def categorical_rgba(data, valid, colors):
    rgba = np.zeros((*data.shape, 4), dtype=np.uint8)

    for code, color in colors.items():
        mask = valid & np.isclose(data, float(code))
        rgba[mask] = hex_rgba(color)

    return rgba


def continuous_rgba(data, valid, minimum, maximum, palette):
    rgba = np.zeros((*data.shape, 4), dtype=np.uint8)
    scaled = np.zeros(data.shape, dtype=np.float32)
    span = maximum - minimum

    if span > 0 and np.any(valid):
        scaled[valid] = np.clip(
            (data[valid] - minimum) / span,
            0.0,
            1.0,
        )

    cmap_name = {
        "ylgn": "YlGn",
        "magma": "magma",
        "viridis": "viridis",
    }.get(str(palette).lower(), "viridis")

    colored = (
        colormaps[cmap_name](scaled) * 255
    ).astype(np.uint8)

    rgba[valid, :3] = colored[valid, :3]
    rgba[valid, 3] = 255

    return rgba


def read_data_and_valid_mask(src):
    """Read band 1 safely regardless of integer or floating-point dtype."""
    band = src.read(1, masked=True)

    # Cast before filling. Integer masked arrays cannot be filled with NaN.
    float_band = band.astype(np.float32)
    data = np.asarray(
        float_band.filled(np.nan),
        dtype=np.float32,
    )

    valid = ~np.ma.getmaskarray(band)
    valid &= np.isfinite(data)

    # Respect GDAL's dataset mask or alpha mask when present.
    dataset_mask = src.dataset_mask()
    valid &= dataset_mask > 0

    # Defensive nodata check for rasters whose mask metadata is incomplete.
    if src.nodata is not None:
        valid &= ~np.isclose(data, float(src.nodata))

    return data, valid


def main():
    if not CATALOG_PATH.exists():
        raise FileNotFoundError(
            f"Catalogue not found: {CATALOG_PATH}"
        )

    PNG_DIR.mkdir(parents=True, exist_ok=True)

    catalog = json.loads(
        CATALOG_PATH.read_text(encoding="utf-8")
    )

    updated_layers = []

    for layer in catalog.get("layers", []):
        source = DOCS / layer["file"]

        if not source.exists():
            print(f"Missing: {source}")
            continue

        with rasterio.open(source) as src:
            data, valid = read_data_and_valid_mask(src)

            if layer["kind"] == "categorical":
                rgba = categorical_rgba(
                    data,
                    valid,
                    catalog["colors"],
                )
            else:
                rgba = continuous_rgba(
                    data,
                    valid,
                    float(layer.get("min", 0.0)),
                    float(layer.get("max", 1.0)),
                    layer.get("palette", "viridis"),
                )

            output_name = f"{source.stem}.png"
            output_path = PNG_DIR / output_name

            Image.fromarray(
                rgba,
                mode="RGBA",
            ).save(
                output_path,
                optimize=True,
            )

            bounds = src.bounds

        updated = dict(layer)
        updated["file"] = f"data_png/{output_name}"
        updated["imageExtent"] = [
            bounds.left,
            bounds.bottom,
            bounds.right,
            bounds.top,
        ]
        updated["projection"] = "EPSG:8857"
        updated_layers.append(updated)

        transparent = int((rgba[:, :, 3] == 0).sum())
        opaque = int((rgba[:, :, 3] == 255).sum())

        print(
            f"Wrote {output_name} | "
            f"transparent={transparent:,} | "
            f"opaque={opaque:,}"
        )

    catalog["layers"] = updated_layers
    catalog["projection"] = "EPSG:8857"
    catalog["rendering"] = "transparent_rgba_png"

    PNG_CATALOG_PATH.write_text(
        json.dumps(
            catalog,
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )

    print(
        f"Wrote {len(updated_layers)} PNG layers and "
        f"{PNG_CATALOG_PATH}"
    )


if __name__ == "__main__":
    main()
