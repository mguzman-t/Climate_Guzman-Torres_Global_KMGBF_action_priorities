#!/usr/bin/env python3
"""
Build transparent Equal Earth PNG layers using high-contrast scientific
colour palettes.

Palette assignment
------------------
Categorical layers
    Use the categorical class colours stored in catalog.json.

Priority-score layers
    Use Plasma.

Winning-score margin
    Use Inferno.

Restoration-count layers
    Use YlGn.

Important
---------
The viewer displays pre-rendered PNG files. Therefore, changing only app.js or
catalogue metadata does not change the map colours. This script must be rerun
whenever the continuous colour palettes are changed.
"""

from pathlib import Path
import json
import shutil

import numpy as np
import rasterio
from PIL import Image
from matplotlib import colormaps


# =============================================================================
# PATHS
# =============================================================================

DOCS = Path("docs")

CATALOG_PATH = (
    DOCS
    / "catalog.json"
)

PNG_DIR = (
    DOCS
    / "data_png"
)

PNG_CATALOG_PATH = (
    DOCS
    / "catalog_png.json"
)


# =============================================================================
# CONFIGURATION
# =============================================================================

# Delete previously rendered PNGs before creating the new versions.
CLEAR_EXISTING_PNGS = True

# Catalogue palette names mapped to Matplotlib registered colour-map names.
PALETTE_NAMES = {
    "plasma": "plasma",
    "inferno": "inferno",
    "ylgn": "YlGn",

    # Retained for backward compatibility.
    "cividis": "cividis",
    "viridis": "viridis",
    "magma": "magma",
}

# Default high-contrast palette for priority-score layers.
DEFAULT_PRIORITY_PALETTE = "plasma"

# Palette used for the difference between winning and runner-up scores.
MARGIN_PALETTE = "inferno"

# Palette used for restoration model-count layers.
RESTORATION_COUNT_PALETTE = "ylgn"


# =============================================================================
# COLOUR UTILITIES
# =============================================================================

def hex_rgba(value):
    """
    Convert a hexadecimal colour to an RGBA tuple with 0-255 channels.

    Examples
    --------
    "#DA7C59"   -> (218, 124, 89, 255)
    "#DA7C59FF" -> (218, 124, 89, 255)
    """

    value = str(value).strip().lstrip("#")

    if len(value) == 6:
        value += "FF"

    if len(value) != 8:
        raise ValueError(
            "Hexadecimal colours must contain 6 or 8 characters. "
            f"Received: {value!r}"
        )

    return tuple(
        int(
            value[index:index + 2],
            16,
        )
        for index in (
            0,
            2,
            4,
            6,
        )
    )


def choose_continuous_palette(layer):
    """
    Select a continuous colour palette from the layer meaning.

    This function deliberately overrides obsolete palette metadata in the
    source catalogue so that rerunning this script always regenerates the
    intended high-contrast PNG colours.
    """

    title = str(
        layer.get(
            "title",
            "",
        )
    ).strip().lower()

    maximum = float(
        layer.get(
            "max",
            1.0,
        )
    )

    # Winning score minus runner-up score.
    if "margin" in title:
        return MARGIN_PALETTE

    # Model count or another multi-level restoration context layer.
    if (
        "restoration" in title
        and maximum > 1.0
    ):
        return RESTORATION_COUNT_PALETTE

    # All normalized and consensus-adjusted action scores.
    return DEFAULT_PRIORITY_PALETTE


def categorical_rgba(
    data,
    valid,
    colors,
):
    """
    Render categorical values using the class colours from catalog.json.

    Pixels outside the valid raster mask remain transparent.
    """

    rgba = np.zeros(
        (
            *data.shape,
            4,
        ),
        dtype=np.uint8,
    )

    for code, color in colors.items():
        numeric_code = float(code)

        mask = (
            valid
            & np.isclose(
                data,
                numeric_code,
            )
        )

        rgba[mask] = hex_rgba(
            color
        )

    return rgba


def continuous_rgba(
    data,
    valid,
    minimum,
    maximum,
    palette,
):
    """
    Render a continuous raster as an opaque RGBA image over valid cells.

    Invalid and NoData pixels remain fully transparent.
    """

    rgba = np.zeros(
        (
            *data.shape,
            4,
        ),
        dtype=np.uint8,
    )

    scaled = np.zeros(
        data.shape,
        dtype=np.float32,
    )

    span = maximum - minimum

    if (
        span > 0
        and np.any(valid)
    ):
        scaled[valid] = np.clip(
            (
                data[valid]
                - minimum
            )
            / span,
            0.0,
            1.0,
        )

    requested = str(
        palette
        or DEFAULT_PRIORITY_PALETTE
    ).lower()

    cmap_name = PALETTE_NAMES.get(
        requested
    )

    if cmap_name is None:
        raise ValueError(
            f"Unsupported palette {palette!r}. "
            f"Supported palettes: "
            f"{sorted(PALETTE_NAMES)}"
        )

    cmap = colormaps.get_cmap(
        cmap_name
    )

    # bytes=True directly returns uint8 channels in the 0-255 range.
    colored = cmap(
        scaled,
        bytes=True,
    ).astype(
        np.uint8
    )

    # Only valid analytical pixels receive visible colours.
    rgba[valid, :3] = colored[valid, :3]

    # Full opacity over valid cells.
    rgba[valid, 3] = 255

    return rgba


# =============================================================================
# RASTER UTILITIES
# =============================================================================

def read_data_and_valid_mask(source):
    """
    Read band 1 as float32 and create a robust valid-data mask.

    The mask combines:

    1. Rasterio's band mask.
    2. Finite-value checking.
    3. The dataset-wide mask.
    4. The explicit NoData value, if defined.
    """

    band = source.read(
        1,
        masked=True,
    )

    # Convert to floating point before filling with NaN.
    float_band = band.astype(
        np.float32
    )

    data = np.asarray(
        float_band.filled(
            np.nan
        ),
        dtype=np.float32,
    )

    valid = (
        ~np.ma.getmaskarray(
            band
        )
    )

    valid &= np.isfinite(
        data
    )

    valid &= (
        source.dataset_mask()
        > 0
    )

    if source.nodata is not None:
        valid &= (
            ~np.isclose(
                data,
                float(
                    source.nodata
                ),
            )
        )

    return (
        data,
        valid,
    )


def clear_existing_outputs():
    """
    Remove previously rendered PNG layers and the old PNG catalogue.
    """

    if not CLEAR_EXISTING_PNGS:
        return

    if PNG_DIR.exists():
        shutil.rmtree(
            PNG_DIR
        )

    if PNG_CATALOG_PATH.exists():
        PNG_CATALOG_PATH.unlink()


# =============================================================================
# LAYER RENDERING
# =============================================================================

def render_layer(
    layer,
    catalog,
):
    """
    Render one catalogue layer and return its updated PNG metadata.
    """

    source_path = (
        DOCS
        / layer["file"]
    )

    if not source_path.exists():
        print(
            f"Missing: {source_path}"
        )

        return None

    with rasterio.open(
        source_path
    ) as source:
        data, valid = (
            read_data_and_valid_mask(
                source
            )
        )

        kind = str(
            layer.get(
                "kind",
                "",
            )
        ).lower()

        if kind == "categorical":
            rgba = categorical_rgba(
                data=data,
                valid=valid,
                colors=catalog["colors"],
            )

            palette_used = "categorical"

        elif kind == "continuous":
            palette_used = (
                choose_continuous_palette(
                    layer
                )
            )

            rgba = continuous_rgba(
                data=data,
                valid=valid,
                minimum=float(
                    layer.get(
                        "min",
                        0.0,
                    )
                ),
                maximum=float(
                    layer.get(
                        "max",
                        1.0,
                    )
                ),
                palette=palette_used,
            )

        else:
            raise ValueError(
                f"Unsupported layer kind "
                f"{kind!r} for "
                f"{source_path}"
            )

        output_name = (
            f"{source_path.stem}.png"
        )

        output_path = (
            PNG_DIR
            / output_name
        )

        image = Image.fromarray(
            rgba,
            mode="RGBA",
        )

        image.save(
            output_path,
            format="PNG",
            optimize=True,
            compress_level=9,
        )

        bounds = source.bounds

    updated = dict(
        layer
    )

    updated["file"] = (
        f"data_png/{output_name}"
    )

    updated["imageExtent"] = [
        float(bounds.left),
        float(bounds.bottom),
        float(bounds.right),
        float(bounds.top),
    ]

    updated["projection"] = (
        "EPSG:8857"
    )

    updated["palette"] = (
        palette_used
    )

    transparent = int(
        (
            rgba[:, :, 3]
            == 0
        ).sum()
    )

    opaque = int(
        (
            rgba[:, :, 3]
            == 255
        ).sum()
    )

    intermediate = int(
        (
            (
                rgba[:, :, 3]
                > 0
            )
            & (
                rgba[:, :, 3]
                < 255
            )
        ).sum()
    )

    print(
        f"Wrote {output_name} | "
        f"palette={palette_used} | "
        f"transparent={transparent:,} | "
        f"opaque={opaque:,} | "
        f"intermediate_alpha={intermediate:,}"
    )

    return updated


# =============================================================================
# MAIN
# =============================================================================

def main():
    if not CATALOG_PATH.exists():
        raise FileNotFoundError(
            CATALOG_PATH
        )

    clear_existing_outputs()

    PNG_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    catalog = json.loads(
        CATALOG_PATH.read_text(
            encoding="utf-8",
        )
    )

    if "layers" not in catalog:
        raise KeyError(
            "catalog.json does not contain "
            "a 'layers' list."
        )

    if "colors" not in catalog:
        raise KeyError(
            "catalog.json does not contain "
            "categorical class colours."
        )

    updated_layers = []

    palette_counts = {}

    for layer in catalog.get(
        "layers",
        [],
    ):
        updated = render_layer(
            layer=layer,
            catalog=catalog,
        )

        if updated is None:
            continue

        updated_layers.append(
            updated
        )

        palette_name = updated.get(
            "palette",
            "unknown",
        )

        palette_counts[
            palette_name
        ] = (
            palette_counts.get(
                palette_name,
                0,
            )
            + 1
        )

    output_catalog = dict(
        catalog
    )

    output_catalog["layers"] = (
        updated_layers
    )

    output_catalog["projection"] = (
        "EPSG:8857"
    )

    output_catalog["rendering"] = (
        "transparent_rgba_png"
    )

    output_catalog["continuous_palettes"] = {
        "priority_scores": (
            DEFAULT_PRIORITY_PALETTE
        ),
        "winning_score_margin": (
            MARGIN_PALETTE
        ),
        "restoration_count": (
            RESTORATION_COUNT_PALETTE
        ),
    }

    PNG_CATALOG_PATH.write_text(
        json.dumps(
            output_catalog,
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )

    print()
    print(
        f"Wrote {len(updated_layers)} "
        f"PNG layers"
    )

    print(
        f"Wrote {PNG_CATALOG_PATH}"
    )

    print(
        "Layer counts by palette:"
    )

    for palette_name in sorted(
        palette_counts
    ):
        print(
            f"  {palette_name}: "
            f"{palette_counts[palette_name]}"
        )


if __name__ == "__main__":
    main()
