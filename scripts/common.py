import numpy as np
import xarray as xr
from config import CRS, FLOAT_NODATA, UINT8_NODATA

def align(data, template):
    return data.reindex_like(template, method='nearest')

def percentile_rank(data):
    values=np.asarray(data.values,dtype=np.float64); valid=np.isfinite(values)
    out=np.full(values.shape,np.nan,dtype=np.float32); x=values[valid]
    if x.size==0: return data.copy(data=out)
    if x.size==1: out[valid]=100.; return data.copy(data=out)
    order=np.argsort(x,kind='mergesort'); sx=x[order]; sr=np.empty(x.size); start=0
    while start<x.size:
        end=start+1
        while end<x.size and sx[end]==sx[start]: end+=1
        sr[start:end]=(start+end-1)/2.; start=end
    ranks=np.empty_like(sr); ranks[order]=sr
    out[valid]=(100*ranks/(x.size-1)).astype(np.float32)
    return data.copy(data=out)

def normalized_percentile(data, mask):
    return (percentile_rank(data.where(mask))/100.).clip(0,1)

def require(dataset,name):
    if name not in dataset:
        raise KeyError(f"Required variable '{name}' is missing from {dataset.encoding.get('source','dataset')}")
    return dataset[name].astype(np.float32)

def write_float(data,path):
    out=data.astype(np.float32).rio.write_crs(CRS).rio.write_nodata(FLOAT_NODATA)
    out.rio.to_raster(path,dtype='float32',compress='LZW',nodata=FLOAT_NODATA)

def write_uint8(data,path):
    out=xr.where(data.notnull(),data,UINT8_NODATA).astype(np.uint8)
    out=out.rio.write_crs(CRS).rio.write_nodata(UINT8_NODATA)
    out.rio.to_raster(path,dtype='uint8',compress='LZW',nodata=UINT8_NODATA)
