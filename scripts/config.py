from pathlib import Path

# Final KM-GBF prioritisation outputs after fixed positive-only pooled
# target-specific calibration. Raw ecological scores and non-winner
# diagnostics are retained alongside the calibrated decision products.
ANALYSIS_ROOT = Path(
    "/capacity/occr_davin/mguzman/P1/Connectivity/approach2/"
    "kmgbf_memberwise_t_calibration_positive_180"
)

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
        "T1: spatial planning under biodiversity and intact ecosystems loss",

    21:
        "T2: ecosystem/connectivity restoration",

    8:
        "T8: climate-adaptation conservation priority",

    10:
        "T10: biodiversity-sensitive production landscapes",

    31:
        "T3.1: protect biodiversity-relevant and climate-persistent diffuse/permeable landscapes",

    32:
        "T3.2: protect concentrated connectivity",

    33:
        "T3.3: support existing PA/OECM under climate vulnerability",

    34:
        "T3.4: safeguard existing PA/OECM under land-use pressure",

    # In the modal map, 97 marks either exact score ties within members
    # or cells where two or more actions share the highest ensemble count.
    97:
        "Tie",
}


ACTION_COLORS = {
    # No positive score
    0:  "#D9D9D9",

    # T1 spatial planning
    1:  "#D55E00",   # vermillion

    # T2 restoration
    21: "#009E73",   # bluish green

    # T8 climate adaptation
    8:  "#CC79A7",   # reddish purple

    # T10 production landscapes
    10: "#E69F00",   # orange

    # T31 diffuse/permeable conservation
    31: "#56B4E9",   # sky blue

    # T32 concentrated connectivity
    32: "#0072B2",   # strong blue

    # T33 existing PA/OECM + climate vulnerability
    33: "#6A3D9A",   # dark purple

    # T34 existing PA/OECM + land-use pressure
    34: "#8C564B",   # brown

    # Tie
    97: "#202020",
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
