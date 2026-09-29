#!/usr/bin/env python3
"""Reproject local Natural Earth countries from EPSG:4326 to EPSG:8857."""
from pathlib import Path
import geopandas as gpd

SOURCE = Path('docs/data_boundaries/ne_110m_admin_0_countries.geojson')
DEST = Path('docs/data_boundaries/ne_110m_admin_0_countries_8857.geojson')

if not SOURCE.exists():
    raise FileNotFoundError(SOURCE)

gdf = gpd.read_file(SOURCE)
if gdf.crs is None:
    gdf = gdf.set_crs('EPSG:4326')
else:
    gdf = gdf.to_crs('EPSG:4326')

projected = gdf.to_crs('EPSG:8857')
DEST.parent.mkdir(parents=True, exist_ok=True)
projected.to_file(DEST, driver='GeoJSON')

print(f'Wrote {DEST}')
print(f'Country features: {len(projected)}')
print(f'Total bounds: {projected.total_bounds.tolist()}')
