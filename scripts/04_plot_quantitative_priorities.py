#!/usr/bin/env python3
"""Create main and action-specific figures from existing score NetCDFs."""
import warnings
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
from config import *
try:
 import cartopy.crs as ccrs; import cartopy.feature as cfeature; HAS_CARTOPY=True
except ImportError: HAS_CARTOPY=False
plt.rcParams.update({'font.size':11,'axes.titlesize':12,'axes.titleweight':'bold','figure.titlesize':18,'savefig.facecolor':'white','figure.facecolor':'white'})

def axes(nr,nc,size):
 return plt.subplots(nr,nc,figsize=size,subplot_kw={'projection':ccrs.PlateCarree()} if HAS_CARTOPY else {},squeeze=False)
def context(ax):
 if HAS_CARTOPY:
  ax.add_feature(cfeature.OCEAN,facecolor='#EAF1F4',zorder=0); ax.add_feature(cfeature.LAND,facecolor='#F3F3F0',zorder=0)
  ax.coastlines(linewidth=.45,color='#454545',zorder=5); ax.add_feature(cfeature.BORDERS,linewidth=.22,edgecolor='#6D6D6D',zorder=5); ax.set_global()
 else: ax.set_xlim(-180,180); ax.set_ylim(-90,90)
def kwargs(): return {'shading':'nearest','rasterized':True,'antialiased':False,'zorder':2,**({'transform':ccrs.PlateCarree()} if HAS_CARTOPY else {})}
def save(fig,name):
 for ext in OUTPUT_FORMATS:
  p=FIGURE_DIR/f'{name}.{ext}'; fig.savefig(p,dpi=DPI,bbox_inches='tight',facecolor='white'); print(f'Wrote {p}')
 plt.close(fig)
def load_mean():
 out={}
 for ssp in SSPS:
  for jump in JUMPS:
   p=MEAN_OUTPUT_DIR/f'kmgbf_action_priority_mean_{ssp}_{jump}.nc'
   if p.exists(): out[(ssp,jump)]=xr.open_dataset(p).load()
 return out

def primary_map(results):
 fig,axs=axes(2,3,FIGSIZE); codes=[*ACTION_CODES,TIE_CODE]; cmap=ListedColormap([ACTION_COLORS[c] for c in codes]); cmap.set_bad((1,1,1,0)); norm=BoundaryNorm(np.arange(-.5,len(codes)+.5,1),cmap.N)
 for r,jump in enumerate(JUMPS):
  for col,ssp in enumerate(SSPS):
   ax=axs[r,col]; ds=results.get((ssp,jump))
   if ds is None: ax.set_axis_off(); continue
   src=ds['primary_action']; enc=xr.full_like(src,np.nan,dtype=np.float32)
   for i,c in enumerate(codes): enc=xr.where(src==c,i,enc)
   ax.pcolormesh(enc.x,enc.y,enc.values,cmap=cmap,norm=norm,**kwargs()); context(ax); ax.set_title(f'{SSP_TITLES[ssp]} | {JUMP_TITLES[jump]}')
 fig.suptitle('Quantitative prioritisation of KM-GBF action classes\nIMAGE-MAgPIE mean',y=.985,fontweight='bold')
 fig.legend(handles=[Patch(facecolor=ACTION_COLORS[c],edgecolor='#333',label=ACTION_LABELS[c]) for c in codes],loc='lower center',bbox_to_anchor=(.5,.008),ncol=5,frameon=False,fontsize=9.5)
 fig.subplots_adjust(left=.015,right=.985,bottom=.16,top=.88,wspace=.02,hspace=.09); save(fig,'Figure_1_quantitative_priority_classes')

def continuous(results,var,title,label,cmap,name,vmax=1):
 fig,axs=axes(2,3,FIGSIZE); mesh=None
 for r,jump in enumerate(JUMPS):
  for col,ssp in enumerate(SSPS):
   ax=axs[r,col]; ds=results.get((ssp,jump))
   if ds is None: ax.set_axis_off(); continue
   d=ds[var]; mesh=ax.pcolormesh(d.x,d.y,d.values,cmap=cmap,vmin=0,vmax=vmax,**kwargs()); context(ax); ax.set_title(f'{SSP_TITLES[ssp]} | {JUMP_TITLES[jump]}')
 fig.suptitle(title,y=.985,fontweight='bold')
 if mesh is not None: cb=fig.colorbar(mesh,ax=axs.ravel().tolist(),orientation='horizontal',fraction=.045,pad=.075,aspect=55); cb.set_label(label,fontsize=12)
 fig.subplots_adjust(left=.015,right=.985,bottom=.12,top=.90,wspace=.02,hspace=.09); save(fig,name)

def atlases():
 for ssp in SSPS:
  for jump in JUMPS:
   p=MEAN_OUTPUT_DIR/f'kmgbf_action_priority_mean_{ssp}_{jump}.nc'
   if not p.exists(): continue
   ds=xr.open_dataset(p).load(); fig,axs=axes(2,4,ATLAS_FIGSIZE); mesh=None
   for ax,code in zip(axs.ravel(),ACTION_CODES):
    d=ds[f'action_{code}_consensus_adjusted_score']; mesh=ax.pcolormesh(d.x,d.y,d.values,cmap='viridis',vmin=0,vmax=1,**kwargs()); context(ax); ax.set_title(ACTION_LABELS[code],fontsize=10)
   fig.suptitle(f'IMAGE-MAgPIE mean priority scores by action\n{SSP_TITLES[ssp]} | {JUMP_TITLES[jump]}',y=.985,fontweight='bold')
   cb=fig.colorbar(mesh,ax=axs.ravel().tolist(),orientation='horizontal',fraction=.04,pad=.07,aspect=60); cb.set_label('Consensus-adjusted within-action percentile score')
   fig.subplots_adjust(left=.015,right=.985,bottom=.12,top=.88,wspace=.03,hspace=.12); save(fig,f'additional_mean_action_scores_{ssp}_{jump}'); ds.close()

def main():
 results=load_mean()
 if not results: raise FileNotFoundError(f'No mean NetCDFs in {MEAN_OUTPUT_DIR}')
 primary_map(results)
 continuous(results,'winning_consensus_adjusted_score','Strength of the quantitatively prioritised action\nIMAGE-MAgPIE mean','Winning consensus-adjusted score','viridis','Figure_2_winning_priority_score')
 continuous(results,'winning_score_margin','Separation between the prioritised action and the alternative\nIMAGE-MAgPIE mean','Winning score minus runner-up score','magma','Figure_3_priority_score_margin')
 continuous(results,'projected_restoration_model_count','Projected restoration context\nNumber of land-use models projecting restoration','Land-use models: 0, 1, or 2','YlGn','Figure_4_projected_restoration_context',2)
 atlases()
 for ds in results.values(): ds.close()
if __name__=='__main__': main()
