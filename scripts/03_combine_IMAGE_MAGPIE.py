#!/usr/bin/env python3
"""Combine IMAGE and MAgPIE action scores and select the quantitative priority."""
import json, warnings
import numpy as np
import xarray as xr
from config import *
from common import align, write_float, write_uint8

def model_path(ssp,lum,jump): return MODEL_OUTPUT_DIR/f'kmgbf_action_scores_{ssp}_{lum}_{jump}.nc'

def choose(scores,valid):
    stack=xr.concat([scores[c] for c in ACTION_CODES],dim='action_code').assign_coords(action_code=ACTION_CODES)
    any_positive=(stack>0).any('action_code'); win=stack.max('action_code').where(any_positive)
    ties=((abs(stack-win)<=SCORE_TIE_TOLERANCE)&(stack>0)).sum('action_code').astype(np.float32)
    vals=stack.fillna(0).values; idx=np.argmax(vals,axis=0)
    primary=valid.copy(data=np.take(np.asarray(ACTION_CODES,dtype=np.uint8),idx)).astype(np.float32)
    primary=xr.where(~any_positive,0,primary); primary=xr.where(ties>1,TIE_CODE,primary).where(valid)
    sorted_vals=np.sort(vals,axis=0); runner=win.copy(data=sorted_vals[-2].astype(np.float32))
    margin=xr.where(any_positive,win-runner,np.nan).astype(np.float32)
    return primary,win,runner,margin,ties

def combine(image,magpie,ssp,jump):
    template=image['valid_data_coverage']; valid=image['valid_data_coverage'].notnull()|align(magpie['valid_data_coverage'],template).notnull()
    out=xr.Dataset(); consensus={}
    for code in ACTION_CODES:
        i=image[f'action_{code}_normalized_score']; m=align(magpie[f'action_{code}_normalized_score'],template)
        ie=image[f'action_{code}_eligible']==1; me=align(magpie[f'action_{code}_eligible'],template)==1
        available=xr.concat([i,m],dim='land_use_model').mean('land_use_model',skipna=True).where(ie|me)
        adjusted=((i.fillna(0)+m.fillna(0))/2).where(ie|me)
        support=ie.astype(np.float32)+me.astype(np.float32)
        out[f'action_{code}_IMAGE_score']=i; out[f'action_{code}_MAgPIE_score']=m
        out[f'action_{code}_mean_available_score']=available
        out[f'action_{code}_consensus_adjusted_score']=adjusted
        out[f'action_{code}_LUM_support_count']=support; consensus[code]=adjusted
    primary,win,runner,margin,ties=choose(consensus,valid)
    out['primary_action']=primary; out['winning_consensus_adjusted_score']=win
    out['runner_up_consensus_adjusted_score']=runner; out['winning_score_margin']=margin; out['winning_tie_count']=ties
    ic=image['projected_restoration_context']; mc=align(magpie['projected_restoration_context'],template)
    out['projected_restoration_context_IMAGE']=ic; out['projected_restoration_context_MAgPIE']=mc
    out['projected_restoration_model_count']=(ic>0).astype(np.float32)+(mc>0).astype(np.float32)
    out['projected_restoration_aligned_with_need_count']=(ic==2).astype(np.float32)+(mc==2).astype(np.float32)
    out.attrs.update({'ssp':ssp,'time_jump':jump,'action_labels':json.dumps(ACTION_LABELS),'primary_selection':'Maximum IMAGE-MAgPIE consensus-adjusted within-action percentile','tie_code':TIE_CODE})
    return out

def save(out,ssp,jump):
    stem=f'{ssp}_IMAGE-MAgPIE-mean_{jump}_native05'; nc=MEAN_OUTPUT_DIR/f'kmgbf_action_priority_mean_{ssp}_{jump}.nc'
    enc={n:{'dtype':'float32','_FillValue':FLOAT_NODATA,'zlib':True,'complevel':4} for n in out.data_vars}; out.to_netcdf(nc,encoding=enc)
    write_uint8(out['primary_action'],MEAN_OUTPUT_DIR/f'primary_action_{stem}.tif')
    write_float(out['winning_consensus_adjusted_score'],MEAN_OUTPUT_DIR/f'winning_score_{stem}.tif')
    write_float(out['winning_score_margin'],MEAN_OUTPUT_DIR/f'winning_score_margin_{stem}.tif')
    write_uint8(out['projected_restoration_model_count'],MEAN_OUTPUT_DIR/f'projected_restoration_model_count_{stem}.tif')
    for code in ACTION_CODES:
      write_float(out[f'action_{code}_mean_available_score'],MEAN_OUTPUT_DIR/f'action_{code}_mean_available_score_{stem}.tif')
      write_float(out[f'action_{code}_consensus_adjusted_score'],MEAN_OUTPUT_DIR/f'action_{code}_consensus_adjusted_score_{stem}.tif')
      write_uint8(out[f'action_{code}_LUM_support_count'],MEAN_OUTPUT_DIR/f'action_{code}_LUM_support_count_{stem}.tif')

def main():
    for ssp in SSPS:
      for jump in JUMPS:
        ip=model_path(ssp,'image',jump); mp=model_path(ssp,'magpie',jump)
        if not ip.exists() or not mp.exists(): warnings.warn(f'Missing IMAGE/MAgPIE files for {ssp} {jump}'); continue
        image=xr.open_dataset(ip).load(); magpie=xr.open_dataset(mp).load(); out=combine(image,magpie,ssp,jump)
        save(out,ssp,jump); image.close(); magpie.close(); out.close(); print(f'Wrote mean {ssp} | {jump}')
if __name__=='__main__': main()
