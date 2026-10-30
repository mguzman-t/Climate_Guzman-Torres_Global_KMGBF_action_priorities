#!/usr/bin/env python3
"""Render PNGs with action-specific agreement colors and corrected restoration scales."""
import json
import shutil
import numpy as np
import rasterio
from PIL import Image
from matplotlib import colormaps
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from config import DOCS, WEB_PNG_DIR

AGREEMENT_CMAP = LinearSegmentedColormap.from_list(
    "ensemble_agreement",
    ["#FFFFFF", "#DEEBF7", "#9ECAE1", "#6BAED6", "#3182BD", "#08519C", "#08306B"]
)


def white_start(base_name):
    base = colormaps.get_cmap(base_name)
    colors = base(np.linspace(0.08, 1.0, 256))
    transition = 28
    target = colors[transition].copy()
    for index in range(transition):
        fraction = index / transition
        colors[index, :3] = (1-fraction) * \
            np.array([1., 1., 1.]) + fraction * target[:3]
    return LinearSegmentedColormap.from_list("white_" + base_name, colors)


def action_cmap(action, catalog):
    color = catalog["colors"][str(action)]
    rgb = to_rgb(color)
    midpoint = tuple(0.62 + 0.38 * channel for channel in rgb)
    return LinearSegmentedColormap.from_list(
        f"action_{action}_agreement", ["#FFFFFF", midpoint, color]
    )


def hexrgba(value):
    value = value.lstrip("#")
    return tuple(int(value[i:i+2], 16) for i in (0, 2, 4)) + (255,)


def choose_cmap(layer, catalog):
    metric = layer.get("metric")
    if metric == "ensemble_agreement" and layer.get("action") is not None:
        return action_cmap(layer["action"], catalog), f"action_{layer['action']}"
    if metric in {"modal_agreement", "connectivity_agreement"}:
        return AGREEMENT_CMAP, "agreement"
    if metric == "restoration_gap":
        return white_start("magma"), "magma"
    if metric in {"recovery_score", "recovery_alignment"}:
        return white_start("YlGn"), "ylgn"
    return white_start("plasma"), "plasma"


def main():
    catalog = json.loads((DOCS / "catalog.json").read_text(encoding="utf-8"))
    if WEB_PNG_DIR.exists():
        shutil.rmtree(WEB_PNG_DIR)
    WEB_PNG_DIR.mkdir(parents=True)
    updated_layers = []

    for layer in catalog["layers"]:
        source_path = DOCS / layer["file"]
        if not source_path.exists():
            print("Missing:", source_path)
            continue
        with rasterio.open(source_path) as source:
            band = source.read(1, masked=True)
            data = band.astype("float32").filled(np.nan)
            valid = np.isfinite(data) & (~np.ma.getmaskarray(band))
            bounds = source.bounds

        rgba = np.zeros((*data.shape, 4), dtype=np.uint8)
        if layer["kind"] == "categorical":
            for code, color in catalog["colors"].items():
                rgba[valid & np.isclose(data, float(code))] = hexrgba(color)
        else:
            minimum = float(layer.get("min", 0))
            maximum = float(layer.get("max", 100))
            visible = valid & (data > 0)
            scaled = np.zeros_like(data, dtype=np.float32)
            if maximum > minimum and np.any(visible):
                scaled[visible] = np.clip(
                    (data[visible] - minimum) / (maximum - minimum), 0, 1)
            cmap, palette = choose_cmap(layer, catalog)
            colors = cmap(scaled, bytes=True).astype(np.uint8)
            rgba[visible, :3] = colors[visible, :3]
            rgba[visible, 3] = 255
            layer = dict(layer)
            layer["palette"] = palette

        output_name = source_path.stem + ".png"
        Image.fromarray(rgba, "RGBA").save(
            WEB_PNG_DIR / output_name, optimize=True, compress_level=9
        )
        updated = dict(layer)
        updated["file"] = f"data_png/{output_name}"
        updated["imageExtent"] = [float(bounds.left), float(
            bounds.bottom), float(bounds.right), float(bounds.top)]
        updated["zeroTransparent"] = layer["kind"] == "continuous"
        updated["valueOpacityEncoding"] = False
        updated_layers.append(updated)
        print("Wrote", output_name)

    catalog["layers"] = updated_layers
    catalog["continuousEncoding"] = "color only; exact zero transparent; all positive values opaque"
    (DOCS / "catalog_png.json").write_text(
        json.dumps(catalog, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"Wrote {len(updated_layers)} PNG layers")


if __name__ == "__main__":
    main()
