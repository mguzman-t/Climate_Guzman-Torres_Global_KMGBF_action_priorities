#!/usr/bin/env python3
"""Convert mean GeoTIFFs to web COGs and write docs/catalog.json."""
import json, shutil, subprocess
from config import *

def convert(src,dst):
 if shutil.which('gdal_translate'):
  subprocess.run(['gdal_translate','-of','COG','-co','COMPRESS=DEFLATE','-co','OVERVIEWS=AUTO',str(src),str(dst)],check=True)
 else:
  shutil.copy2(src,dst)

def main():
 layers=[]
 for ssp in SSPS:
  for jump in JUMPS:
   stem=f'{ssp}_IMAGE-MAgPIE-mean_{jump}_native05'
   specs=[(f'primary_action_{stem}.tif',f'primary_{ssp}_{jump}.tif','Prioritised action','categorical',{}),
          (f'winning_score_{stem}.tif',f'winning_{ssp}_{jump}.tif','Winning priority score','continuous',{'min':0,'max':1,'palette':'viridis'}),
          (f'winning_score_margin_{stem}.tif',f'margin_{ssp}_{jump}.tif','Winning score margin','continuous',{'min':0,'max':1,'palette':'magma'}),
          (f'projected_restoration_model_count_{stem}.tif',f'restoration_{ssp}_{jump}.tif','Projected restoration model count','continuous',{'min':0,'max':2,'palette':'ylgn'})]
   for code in ACTION_CODES:
    specs.append((f'action_{code}_consensus_adjusted_score_{stem}.tif',f'action_{code}_{ssp}_{jump}.tif',ACTION_LABELS[code],'continuous',{'min':0,'max':1,'palette':'viridis','action':code}))
   for srcname,outname,title,kind,extra in specs:
    src=MEAN_OUTPUT_DIR/srcname
    if not src.exists(): print(f'Missing {src}'); continue
    convert(src,WEB_DATA_DIR/outname); layers.append({'file':f'data/{outname}','title':title,'kind':kind,'ssp':ssp,'period':jump,**extra})
 catalog={'classes':{str(k):v for k,v in ACTION_LABELS.items() if k!=0},'colors':{str(k):v for k,v in ACTION_COLORS.items()},'layers':layers}
 (WEB_DATA_DIR.parent/'catalog.json').write_text(json.dumps(catalog,indent=2),encoding='utf-8')
 print(f'Wrote {len(layers)} web layers')
if __name__=='__main__': main()
