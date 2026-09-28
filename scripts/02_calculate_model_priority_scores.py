#!/usr/bin/env python3
"""Calculate eligibility, raw scores, and within-action percentiles for each LUM."""
import json, warnings
import numpy as np
import xarray as xr
from config import *
from common import normalized_percentile, require, write_float, write_uint8

def input_path(ssp,lum,jump):
    return INPUT_DIR/f'kmgbf_score_selected_inputs_{ssp}_{lum}_{jump}.nc'

def calculate(ds):
    valid=require(ds,'valid_data_coverage').notnull()
    high=require(ds,'high_biodiversity_value')==1; broad=require(ds,'broad_biodiversity_value')==1
    stable=require(ds,'stable_for_target3')==1; declining8=require(ds,'declining_for_target8')==1
    declining=require(ds,'declining_richness')==1
    mean_decline=require(ds,'mean_decline_for_existing_coverage')==1
    c=require(ds,'connectivity_class'); lu=require(ds,'land_use_transition')
    coverage=require(ds,'pa_oecm_fraction').clip(0,1)
    gap=(coverage<PA_OECM_GAP_THRESHOLD)&valid; covered=(coverage>=PA_OECM_GAP_THRESHOLD)&valid
    imp=(c==1)&valid; dif=(c==2)&valid; inte=(c==3)&valid; cha=(c==4)&valid
    lim=(c==5)&valid; per=(c==6)&valid; deg=(c==7)&valid
    need=imp|lim|deg; regional=inte|cha
    stable_lu=(lu==5)&valid; projected_restoration=(lu==4)&valid
    prod=lu.isin([1,2,3])&valid; conv=lu.isin([0,1,2,3])&valid
    eligibility={1:high&conv&valid,2:high&need&valid,8:high&declining8&stable_lu&valid,
                 10:broad&prod&valid,31:gap&high&stable&dif&valid,
                 32:gap&high&stable&regional&valid,34:gap&broad&stable&per&valid,
                 33:covered&high&declining&mean_decline&conv&valid}
    biod=(require(ds,'vertebrate_value_index')/100).clip(0,1)
    breadth=require(ds,'declining_taxon_fraction').clip(0,1)
    mean_change=require(ds,'mean_taxon_richness_change_percent')
    mean_decline_score=normalized_percentile((-mean_change).clip(min=0),valid)
    mean_persistence=normalized_percentile(mean_change,valid)
    loss=(breadth+mean_decline_score)/2
    persistence=((1-breadth).clip(0,1)+mean_persistence)/2
    gap_score=((PA_OECM_GAP_THRESHOLD-coverage)/PA_OECM_GAP_THRESHOLD).clip(0,1)
    need_score=xr.where(deg,1,xr.where(imp,.85,xr.where(lim,.70,0))).astype(np.float32)
    conn_score=xr.where(cha,1,xr.where(inte,.85,xr.where(per,.70,xr.where(dif,.60,0)))).astype(np.float32)
    conv_score=normalized_percentile(require(ds,'conversion_magnitude'),conv)
    prod_score=normalized_percentile(require(ds,'production_magnitude'),prod)
    raw={
      1:WEIGHTS[1]['biodiversity']*biod+WEIGHTS[1]['conversion']*conv_score+WEIGHTS[1]['coverage_gap']*gap_score+WEIGHTS[1]['richness_loss']*loss,
      2:WEIGHTS[2]['biodiversity']*biod+WEIGHTS[2]['connectivity_need']*need_score+WEIGHTS[2]['richness_loss']*loss,
      8:WEIGHTS[8]['biodiversity']*biod+WEIGHTS[8]['decline_breadth']*breadth+WEIGHTS[8]['mean_decline']*mean_decline_score,
      10:WEIGHTS[10]['biodiversity']*biod+WEIGHTS[10]['production_change']*prod_score+WEIGHTS[10]['richness_loss']*loss,
      31:WEIGHTS[31]['biodiversity']*biod+WEIGHTS[31]['persistence']*persistence+WEIGHTS[31]['coverage_gap']*gap_score,
      32:WEIGHTS[32]['biodiversity']*biod+WEIGHTS[32]['connectivity']*conn_score+WEIGHTS[32]['persistence']*persistence+WEIGHTS[32]['coverage_gap']*gap_score,
      34:WEIGHTS[34]['biodiversity']*biod+WEIGHTS[34]['connectivity']*conn_score+WEIGHTS[34]['persistence']*persistence+WEIGHTS[34]['coverage_gap']*gap_score,
      33:WEIGHTS[33]['biodiversity']*biod+WEIGHTS[33]['richness_loss']*loss+WEIGHTS[33]['conversion']*conv_score,
    }
    out=xr.Dataset(); out['valid_data_coverage']=xr.where(valid,1,np.nan)
    context=xr.zeros_like(lu,dtype=np.float32)
    context=xr.where(projected_restoration&~need,1,context)
    context=xr.where(projected_restoration&need,2,context).where(valid)
    out['projected_restoration_context']=context
    out['connectivity_restoration_need']=need.astype(np.float32).where(valid)
    out['richness_loss_score']=loss.where(valid); out['coverage_gap_score']=gap_score.where(valid)
    out['conversion_magnitude_score']=conv_score.where(valid)
    for code in ACTION_CODES:
        out[f'action_{code}_eligible']=eligibility[code].astype(np.float32).where(valid)
        out[f'action_{code}_raw_score']=raw[code].where(eligibility[code])
        out[f'action_{code}_normalized_score']=normalized_percentile(raw[code],eligibility[code]).where(eligibility[code])
    return out

def save(out,ssp,lum,jump):
    stem=f'{ssp}_{lum}_{jump}_native05'; nc=MODEL_OUTPUT_DIR/f'kmgbf_action_scores_{ssp}_{lum}_{jump}.nc'
    enc={n:{'dtype':'float32','_FillValue':FLOAT_NODATA,'zlib':True,'complevel':4} for n in out.data_vars}
    out.to_netcdf(nc,encoding=enc)
    for code in ACTION_CODES:
        write_float(out[f'action_{code}_normalized_score'],MODEL_OUTPUT_DIR/f'action_{code}_normalized_score_{stem}.tif')
        write_uint8(out[f'action_{code}_eligible'],MODEL_OUTPUT_DIR/f'action_{code}_eligible_{stem}.tif')
    write_uint8(out['projected_restoration_context'],MODEL_OUTPUT_DIR/f'projected_restoration_context_{stem}.tif')

def main():
    for ssp in SSPS:
      for jump in JUMPS:
       for lum in LUMS:
        path=input_path(ssp,lum,jump)
        if not path.exists(): warnings.warn(f'Missing {path}'); continue
        ds=xr.open_dataset(path).load(); out=calculate(ds)
        out.attrs.update({'ssp':ssp,'land_use_model':lum,'time_jump':jump,'action_labels':json.dumps(ACTION_LABELS),'weights':json.dumps(WEIGHTS)})
        save(out,ssp,lum,jump); ds.close(); out.close(); print(f'Wrote {ssp} | {lum} | {jump}')
if __name__=='__main__': main()
