#!/usr/bin/env python3
"""
Render web PNG layers from catalog.json.

Continuous layers use colour only:
- exact zero is transparent;
- all positive finite values are opaque;
- values are scaled using the limits recorded in catalog.json.

Action-selection agreement uses the colour assigned to each KM-GBF action.
Raw ecological scores and calibrated scores use distinct continuous palettes
so that their different scientific meanings remain visually explicit.
"""

from __future__ import annotations

import json
import shutil

import numpy as np
import rasterio

from PIL import Image
from matplotlib import colormaps
from matplotlib.colors import LinearSegmentedColormap, to_rgb

from config import (
    DOCS,
    WEB_PNG_DIR,
)


AGREEMENT_CMAP = LinearSegmentedColormap.from_list(
    "ensemble_agreement",
    [
        "#FFFFFF",
        "#DEEBF7",
        "#9ECAE1",
        "#6BAED6",
        "#3182BD",
        "#08519C",
        "#08306B",
    ],
)


def white_start(base_name):
    """
    Return a continuous colour map beginning close to white.

    This keeps low positive values visually distinct from transparent zero.
    """
    base = colormaps.get_cmap(base_name)

    colors = base(
        np.linspace(
            0.08,
            1.0,
            256,
        )
    )

    transition = 28
    target = colors[transition].copy()

    for index in range(transition):
        fraction = index / transition

        colors[index, :3] = (
            (1.0 - fraction)
            * np.array(
                [
                    1.0,
                    1.0,
                    1.0,
                ]
            )
            + fraction
            * target[:3]
        )

    return LinearSegmentedColormap.from_list(
        "white_" + base_name,
        colors,
    )


def action_cmap(action, catalog):
    """Build a white-to-action-colour map for selection agreement."""
    color = catalog["colors"][str(action)]

    rgb = to_rgb(color)

    midpoint = tuple(
        0.62 + 0.38 * channel
        for channel in rgb
    )

    return LinearSegmentedColormap.from_list(
        f"action_{action}_agreement",
        [
            "#FFFFFF",
            midpoint,
            color,
        ],
    )


def hexrgba(value):
    """Convert a hexadecimal RGB colour to an opaque RGBA tuple."""
    value = value.lstrip("#")

    red = int(
        value[0:2],
        16,
    )

    green = int(
        value[2:4],
        16,
    )

    blue = int(
        value[4:6],
        16,
    )

    return (
        red,
        green,
        blue,
        255,
    )


def choose_cmap(layer, catalog):
    """Choose the colour map appropriate for each continuous metric."""
    metric = layer.get("metric")

    if (
        metric == "ensemble_agreement"
        and layer.get("action") is not None
    ):
        action = layer["action"]

        return (
            action_cmap(
                action,
                catalog,
            ),
            f"action_{action}",
        )

    if metric in {
        "modal_agreement",
        "connectivity_agreement",
    }:
        return (
            AGREEMENT_CMAP,
            "agreement",
        )

    if metric == "mean_calibrated_score":
        return (
            white_start("plasma"),
            "plasma",
        )

    if metric == "mean_raw_common_score":
        return (
            white_start("viridis"),
            "viridis",
        )

    if metric == "winning_margin":
        return (
            white_start("plasma"),
            "plasma",
        )

    if metric == "available_member_count":
        return (
            white_start("viridis"),
            "viridis",
        )

    if metric == "restoration_gap":
        return (
            white_start("magma"),
            "magma",
        )

    if metric in {
        "recovery_score",
        "recovery_alignment",
    }:
        return (
            white_start("YlGn"),
            "ylgn",
        )

    return (
        white_start("plasma"),
        "plasma",
    )


def main():
    """Render all catalogue layers and create catalog_png.json."""
    catalog_path = (
        DOCS
        / "catalog.json"
    )

    if not catalog_path.is_file():
        raise FileNotFoundError(
            catalog_path
        )

    catalog = json.loads(
        catalog_path.read_text(
            encoding="utf-8",
        )
    )

    if WEB_PNG_DIR.exists():
        shutil.rmtree(
            WEB_PNG_DIR
        )

    WEB_PNG_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    updated_layers = []

    for layer in catalog["layers"]:
        source_path = (
            DOCS
            / layer["file"]
        )

        if not source_path.is_file():
            print(
                "Missing:",
                source_path,
            )
            continue

        with rasterio.open(source_path) as source:
            band = source.read(
                1,
                masked=True,
            )

            data = (
                band
                .astype("float32")
                .filled(np.nan)
            )

            valid = (
                np.isfinite(data)
                & (
                    ~np.ma.getmaskarray(
                        band
                    )
                )
            )

            bounds = source.bounds

        rgba = np.zeros(
            (
                data.shape[0],
                data.shape[1],
                4,
            ),
            dtype=np.uint8,
        )

        if layer["kind"] == "categorical":

            for code, color in catalog["colors"].items():
                class_mask = (
                    valid
                    & np.isclose(
                        data,
                        float(code),
                    )
                )

                rgba[class_mask] = hexrgba(
                    color
                )

        else:
            minimum = float(
                layer.get(
                    "min",
                    0,
                )
            )

            maximum = float(
                layer.get(
                    "max",
                    100,
                )
            )

            visible = (
                valid
                & (data > 0)
            )

            scaled = np.zeros_like(
                data,
                dtype=np.float32,
            )

            if (
                maximum > minimum
                and np.any(visible)
            ):
                scaled[visible] = np.clip(
                    (
                        data[visible]
                        - minimum
                    )
                    / (
                        maximum
                        - minimum
                    ),
                    0.0,
                    1.0,
                )

            cmap, palette = choose_cmap(
                layer,
                catalog,
            )

            colors = cmap(
                scaled,
                bytes=True,
            ).astype(
                np.uint8
            )

            rgba[
                visible,
                :3,
            ] = colors[
                visible,
                :3,
            ]

            rgba[
                visible,
                3,
            ] = 255

            layer = dict(
                layer
            )

            layer["palette"] = palette

        output_name = (
            source_path.stem
            + ".png"
        )

        output_path = (
            WEB_PNG_DIR
            / output_name
        )

        Image.fromarray(
            rgba
        ).save(
            output_path,
            optimize=True,
            compress_level=9,
        )

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

        updated["zeroTransparent"] = (
            layer["kind"]
            == "continuous"
        )

        updated["valueOpacityEncoding"] = False

        updated_layers.append(
            updated
        )

        print(
            "Wrote",
            output_name,
        )

    catalog["layers"] = updated_layers

    catalog["continuousEncoding"] = (
        "color only; exact zero transparent; "
        "all positive values opaque"
    )

    output_catalog = (
        DOCS
        / "catalog_png.json"
    )

    output_catalog.write_text(
        json.dumps(
            catalog,
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )

    print()
    print(
        f"Wrote {len(updated_layers)} PNG layers"
    )

    print(
        f"PNG catalogue: {output_catalog}"
    )


if __name__ == "__main__":
    main()