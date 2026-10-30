# KM-GBF continuous member-wise viewer

Viewer/export layer for the 50-km Equal Earth KM-GBF prioritisation outputs.

## Pipeline

```bash
cd scripts
python validate_outputs.py
python build_web_catalog.py
python render_png_layers.py
```

Serve `docs/` locally before publishing:

```bash
python3 -m http.server 8000 --directory docs
```

Then open `http://localhost:8000`.

The analytical results are not recomputed here. This repository consumes the finalized SSP x radius outputs from `kmgbf_memberwise_continuous`.
