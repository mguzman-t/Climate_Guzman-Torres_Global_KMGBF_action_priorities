#!/usr/bin/env python3
from pathlib import Path
from urllib.request import Request, urlopen
import json

URL = 'https://cdn.jsdelivr.net/gh/nvkelso/natural-earth-vector@master/geojson/ne_110m_admin_0_countries.geojson'
DEST = Path('docs/data_boundaries/ne_110m_admin_0_countries.geojson')
DEST.parent.mkdir(parents=True, exist_ok=True)
request = Request(URL, headers={'User-Agent': 'KMGBF-Equal-Earth-viewer/1.0'})
with urlopen(request, timeout=60) as response:
    content = response.read()
data = json.loads(content.decode('utf-8'))
features = data.get('features', [])
if not features:
    raise RuntimeError('Downloaded boundary file contains no features')
DEST.write_bytes(content)
print(f'Wrote {DEST}')
print(f'Country features: {len(features)}')
