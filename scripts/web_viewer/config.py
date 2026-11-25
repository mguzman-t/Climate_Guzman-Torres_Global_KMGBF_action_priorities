import os
from pathlib import Path


# Final KM-GBF prioritisation outputs after fixed positive-only pooled
# target-specific calibration.
#
# KMGBF_ANALYSIS_ROOT can be set externally when the analysis outputs are
# stored elsewhere. The Ubelix path is retained as an example.
ANALYSIS_ROOT = Path(
    os.environ.get(
        "KMGBF_ANALYSIS_ROOT",
        (
            "/capacity/occr_davin/mguzman/P1/Connectivity/approach2/"
            "kmgbf_memberwise_t_calibration_positive_180"
        ),
    )
).expanduser().resolve()

# Expected ensemble size per SSP x radius:
# 2 connectivity taxa x 3 SDMs x 5 GCMs x 2 LUMs x 3 resistance curves
EXPECTED_MEMBER_COUNT = 180

REPO_ROOT = Path(__file__).resolve().parents[1]

DOCS = REPO_ROOT / "docs"
WEB_DATA_DIR = DOCS / "data"
WEB_PNG_DIR = DOCS / "data_png"
FIGURE_DIR = REPO_ROOT / "figures"


SSPS = [
    "ssp126",
    "ssp370",
    "ssp585",
]

RADII = [
    6,
    14,
    20,
]


SSP_LABELS = {
    "ssp126": "SSP1-2.6",
    "ssp370": "SSP3-7.0",
    "ssp585": "SSP5-8.5",
}

RADIUS_LABELS = {
    6: "300 km",
    14: "700 km",
    20: "1,000 km",
}


ACTION_CODES = [
    1,
    21,
    8,
    10,
    31,
    32,
    33,
    34,
]

TIE_CODE = 97


ACTION_LABELS = {
    0:
        "No positive policy score",

    1:
        "T1: Spatial planning in vulnerable landscapes",

    21:
        "T2: Restore vertebrates-relevant degraded landscapes",

    8:
        "T8: Promote climate-adaptation",

    10:
        "T10: Enhance sustainability in agricultural landscapes",

    31:
        "T3.1: Protect diverse, well-connected and climate-resilient landscapes",

    32:
        "T3.2: Protect landscapes important for connectivity",

    33:
        "T3.3: Support existing PA/OECM under climate vulnerability",

    34:
        "T3.4: Safeguard existing PA/OECM under land-use pressure",

    # In the modal map, 97 marks either exact score ties within members
    # or cells where two or more actions share the highest ensemble count.
    97:
        "Tie",
}


ACTION_COLORS = {
    0: "#B8B8B8",
    1: "#D55E00",
    21: "#009E73",
    8: "#CC79A7",
    10: "#E5C000",
    31: "#56B4E9",
    32: "#006D77",
    33: "#6F4C9B",
    34: "#8C564B",
    97: "#4D4D4D",
}


for path in [
    WEB_DATA_DIR,
    WEB_PNG_DIR,
    FIGURE_DIR,
]:
    path.mkdir(
        parents=True,
        exist_ok=True,
    )
