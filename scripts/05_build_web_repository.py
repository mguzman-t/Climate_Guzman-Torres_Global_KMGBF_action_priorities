#!/usr/bin/env python3
"""Build EPSG:8857 Equal Earth COGs and web catalog."""
import json, shutil, subprocess
from pathlib import Path
from config import ACTION_CODES, ACTION_COLORS, ACTION_LABELS, JUMPS, LUMS, MEAN_OUTPUT_DIR, MODEL_OUTPUT_DIR, SSPS, WEB_DATA_DIR

def convert(source: Path, destination: Path, kind: str):
    destination.parent.mkdir(parents=True, exist_ok=True)
    if not shutil.which('gdalwarp') or not shutil.which('gdal_translate'):
        raise RuntimeError('gdalwarp and gdal_translate are required')
    tmp=destination.with_suffix('.tmp.tif')
    resampling='near' if kind=='categorical' else 'bilinear'
    nodata='255' if kind=='categorical' else '-9999'
    subprocess.run(['gdalwarp','-overwrite','-t_srs','EPSG:8857','-r',resampling,'-srcnodata',nodata,'-dstnodata',nodata,'-multi','-wo','NUM_THREADS=ALL_CPUS','-co','TILED=YES','-co','COMPRESS=DEFLATE',str(source),str(tmp)],check=True)
    subprocess.run(['gdal_translate','-of','COG','-co','COMPRESS=DEFLATE','-co','OVERVIEWS=AUTO',str(tmp),str(destination)],check=True)
    tmp.unlink(missing_ok=True)

def add(layers, source, filename, title, kind, ssp, period, model, **extra):
    if not source.exists():
        print(f'Missing: {source}')
        return
    convert(source, WEB_DATA_DIR/filename, kind)
    layers.append({'file':f'data/{filename}','title':title,'kind':kind,'ssp':ssp,'period':period,'model':model,**extra})

def main():
    WEB_DATA_DIR.mkdir(parents=True,exist_ok=True)
    layers=[]
    for ssp in SSPS:
      for jump in JUMPS:
        mean=f'{ssp}_IMAGE-MAgPIE-mean_{jump}_native05'
        add(layers,MEAN_OUTPUT_DIR/f'primary_action_{mean}.tif',f'primary_{ssp}_{jump}.tif','Prioritised action','categorical',ssp,jump,'combined')
        add(layers,MEAN_OUTPUT_DIR/f'winning_score_{mean}.tif',f'winning_{ssp}_{jump}.tif','Winning priority score','continuous',ssp,jump,'combined',min=0,max=1,palette='viridis')
        add(layers,MEAN_OUTPUT_DIR/f'winning_score_margin_{mean}.tif',f'margin_{ssp}_{jump}.tif','Winning score margin','continuous',ssp,jump,'combined',min=0,max=1,palette='magma')
        for code in ACTION_CODES:
          add(layers,MEAN_OUTPUT_DIR/f'action_{code}_consensus_adjusted_score_{mean}.tif',f'action_{code}_combined_{ssp}_{jump}.tif',ACTION_LABELS[code],'continuous',ssp,jump,'combined',action=code,min=0,max=1,palette='viridis')
        for lum in LUMS:
          stem=f'{ssp}_{lum}_{jump}_native05'
          for code in ACTION_CODES:
            add(layers,MODEL_OUTPUT_DIR/f'action_{code}_normalized_score_{stem}.tif',f'action_{code}_{lum}_{ssp}_{jump}.tif',ACTION_LABELS[code],'continuous',ssp,jump,lum,action=code,min=0,max=1,palette='viridis')
    catalog={'classes':{str(k):v for k,v in ACTION_LABELS.items() if k!=0},'colors':{str(k):v for k,v in ACTION_COLORS.items()},'model_labels':{'combined':'IMAGE–MAgPIE combined','image':'IMAGE','magpie':'MAgPIE'},'projection':'EPSG:8857','layers':layers}
    (WEB_DATA_DIR.parent/'catalog.json').write_text(json.dumps(catalog,indent=2,ensure_ascii=False),encoding='utf-8')
    print(f'Wrote {len(layers)} Equal Earth web layers')
if __name__=='__main__': main()
