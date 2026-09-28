from pathlib import Path

INPUT_DIR = Path('/capacity/occr_davin/mguzman/P1/Connectivity/kmgbf_05grid_quant_score')
OUTPUT_DIR = Path('/capacity/occr_davin/mguzman/P1/Connectivity/kmgbf_05grid_quant_priority')
MODEL_OUTPUT_DIR = OUTPUT_DIR / 'model_scores'
MEAN_OUTPUT_DIR = OUTPUT_DIR / 'IMAGE_MAGPIE_mean'
FIGURE_DIR = OUTPUT_DIR / 'figures'
WEB_DATA_DIR = Path(__file__).resolve().parents[1] / 'docs' / 'data'

for directory in [OUTPUT_DIR, MODEL_OUTPUT_DIR, MEAN_OUTPUT_DIR, FIGURE_DIR, WEB_DATA_DIR]:
    directory.mkdir(parents=True, exist_ok=True)

SSPS = ['ssp126soc-adapt', 'ssp370soc-adapt', 'ssp585soc-adapt']
LUMS = ['image', 'magpie']
JUMPS = ['J1', 'J2']
SSP_TITLES = {'ssp126soc-adapt':'SSP1-2.6','ssp370soc-adapt':'SSP3-7.0','ssp585soc-adapt':'SSP5-8.5'}
LUM_TITLES = {'image':'IMAGE','magpie':'MAgPIE'}
JUMP_TITLES = {'J1':'2041-2070','J2':'2071-2100'}

ACTION_CODES = [1, 2, 8, 10, 31, 32, 34, 33]
TIE_CODE = 97
SCORE_TIE_TOLERANCE = 1e-6
PA_OECM_GAP_THRESHOLD = 0.30
ACTION_LABELS = {
    0:'No eligible action', 1:'T1: spatial planning',
    2:'T2: connectivity restoration',
    8:'T8: climate vulnerability',
    10:'T10: biodiversity-sensitive production landscapes',
    31:'T3: climate-resilient conservation',
    32:'T3: ecological OECM screening',
    34:'T3: local area-based conservation',
    33:'T3: avoid encroachment', 97:'Unresolved score tie',
}
ACTION_COLORS = {1:'#DA7C59',2:'#3D8B5E',8:'#BE7EA9',10:'#EEBA4D',31:'#87CEEB',32:'#2A75AD',34:'#30B1B5',33:'#4A3991',97:'#202020'}
WEIGHTS = {
  1:{'biodiversity':.25,'conversion':.30,'coverage_gap':.20,'richness_loss':.25},
  2:{'biodiversity':1/3,'connectivity_need':1/3,'richness_loss':1/3},
  8:{'biodiversity':1/3,'decline_breadth':1/3,'mean_decline':1/3},
  10:{'biodiversity':1/3,'production_change':1/3,'richness_loss':1/3},
  31:{'biodiversity':1/3,'persistence':1/3,'coverage_gap':1/3},
  32:{'biodiversity':.25,'connectivity':.25,'persistence':.25,'coverage_gap':.25},
  34:{'biodiversity':.25,'connectivity':.25,'persistence':.25,'coverage_gap':.25},
  33:{'biodiversity':1/3,'richness_loss':1/3,'conversion':1/3},
}
CRS='EPSG:4326'; FLOAT_NODATA=-9999.0; UINT8_NODATA=255
DPI=300; FIGSIZE=(20,10.2); ATLAS_FIGSIZE=(18,12.5); OUTPUT_FORMATS=['png','pdf']
