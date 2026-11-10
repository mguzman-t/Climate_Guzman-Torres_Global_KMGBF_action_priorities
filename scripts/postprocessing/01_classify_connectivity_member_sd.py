#!/usr/bin/env python3
"""
Recreate and export the 180 taxon-specific member-level connectivity classes.

This script keeps taxon as an ensemble dimension.
Each taxon-specific flow surface is classified with its own spatial mean/SD:

    connectivity_class[member, row, col]

    member = taxon x sdm x gcm x lum x curve = 2 x 3 x 5 x 2 x 3 = 180

Class codes:
    1 impeded
    2 diffuse
    3 channelized
    4 intensified
    255 NoData

Output location:
    connectivity_agreement_sd4_ee50km/<ssp>_rad<radius>/
    connectivity_classes_by_taxon_<ssp>_rad<radius>.npz

Member identifiers:
    member_id        taxon__sdm__gcm__lum__curve
    pooled_member_id sdm__gcm__lum__curve   (matches the 90-member archive)

The common valid footprint requires finite values in all 180 surfaces.
"""

import argparse
import itertools
import json
import re
from pathlib import Path

import numpy as np
import rasterio

BASE = Path('/capacity/occr_davin/mguzman/P1/Connectivity/approach2')
INPUT_ROOT = BASE / 'results_omniscape_ee50km_Belote'
OUTPUT_ROOT = BASE / 'connectivity_agreement_sd4_ee50km/'
TEMPLATE = (
    BASE / 'gross_change_ee50km/bird/members/J1/ssp126'
    / 'bird_gam_gfdl-esm4_ssp126_J1_L_ee50km.tif'
)

SSPS = ['ssp126', 'ssp370', 'ssp585']
RADII = [6, 14, 20]
TAXA = ['bird', 'mammal']
SDMS = ['gam', 'glm', 'randomforest']
GCMS = [
    'gfdl-esm4', 'ipsl-cm6a-lr', 'mpi-esm1-2-hr',
    'mri-esm2-0', 'ukesm1-0-ll',
]
LUMS = ['image', 'magpie']
CURVES = {'sensitive': 'neg4', 'linear': '0p25', 'tolerant': '4'}
NODATA_CLASS = 255

PATTERN = re.compile(
    r'^(?P<taxon>bird|mammal)_'
    r'(?P<sdm>gam|glm|randomforest)_'
    r'(?P<gcm>gfdl-esm4|ipsl-cm6a-lr|mpi-esm1-2-hr|mri-esm2-0|ukesm1-0-ll)_'
    r'(?P<ssp>ssp126|ssp370|ssp585)_'
    r'(?P<lum>image|magpie)_'
    r'(?P<curve>sensitive|linear|tolerant)_c(?P<c>neg4|0p25|4)_'
    r'(?:(?P<power>hpow[^_]+)_)?waterR50_'
    r'J1_rad(?P<radius>6|14|20)(?:_RpresStrength)?$'
)


def load_template():
    with rasterio.open(TEMPLATE) as src:
        values = src.read(1, masked=True).astype(float).filled(np.nan)
        target = {
            'crs': src.crs,
            'transform': src.transform,
            'shape': src.shape,
            'footprint': np.isfinite(values),
        }
    if target['crs'] != rasterio.crs.CRS.from_epsg(8857):
        raise ValueError('Template must use EPSG:8857')
    return target


def read_raster(path, target):
    with rasterio.open(path) as src:
        if (
            src.count != 1
            or src.crs != target['crs']
            or src.transform != target['transform']
            or src.shape != target['shape']
        ):
            raise ValueError(f'Grid mismatch: {path}')
        array = src.read(1, masked=True).astype(float).filled(np.nan)
    finite = np.isfinite(array)
    if not finite.any():
        raise ValueError(f'No finite values: {path}')
    if np.any(array[finite] < 0):
        raise ValueError(f'Negative normalized flow: {path}')
    return array


def discover():
    index = {}
    for path in sorted(INPUT_ROOT.rglob('normalized_cum_currmap.tif')):
        match = PATTERN.fullmatch(path.parent.name)
        if match is None:
            continue
        d = match.groupdict()
        if d['c'] != CURVES[d['curve']]:
            continue
        key = (
            d['taxon'], d['sdm'], d['gcm'], d['ssp'], d['lum'],
            d['curve'], int(d['radius']),
        )
        if key in index:
            raise ValueError(f'Duplicate connectivity run: {key}')
        index[key] = path
    if not index:
        raise FileNotFoundError(f'No matching runs under {INPUT_ROOT}')
    return index


def members():
    # Order: taxon outermost, so members 0-89 are bird and 90-179 are mammal,
    # each block following the same sdm/gcm/lum/curve order as the pooled run.
    return list(itertools.product(TAXA, SDMS, GCMS, LUMS, CURVES))


def member_id(member):
    return '__'.join(member)


def pooled_member_id(member):
    return '__'.join(member[1:])


def check_completeness(index, ssp, radius):
    missing = []
    for taxon, sdm, gcm, lum, curve in members():
        key = (taxon, sdm, gcm, ssp, lum, curve, radius)
        if key not in index:
            missing.append(key)
    if missing:
        listing = '\n'.join(str(k) for k in missing)
        raise FileNotFoundError(
            f'{len(missing)} missing connectivity runs for {ssp} rad{radius}:\n'
            f'{listing}'
        )


def member_surface(index, member, ssp, radius, target):
    taxon, sdm, gcm, lum, curve = member
    key = (taxon, sdm, gcm, ssp, lum, curve, radius)
    if key not in index:
        raise FileNotFoundError(f'Missing connectivity run: {key}')
    flow = read_raster(index[key], target)
    flow[~target['footprint']] = np.nan
    return flow


def classify_member(flow, common):
    valid = common & np.isfinite(flow)
    values = flow[valid]
    if values.size < 2:
        raise ValueError('Member has fewer than two valid spatial cells')

    mu = float(np.mean(values))
    sigma = float(np.std(values, ddof=0))
    if not np.isfinite(sigma) or sigma <= 0:
        raise ValueError('Member spatial SD is zero or invalid')

    impeded_max = mu - 0.5 * sigma
    diffuse_max = mu + sigma
    channelized_max = mu + 2.0 * sigma

    classes = np.full(flow.shape, NODATA_CLASS, dtype=np.uint8)
    classes[valid & (flow < impeded_max)] = 1
    classes[valid & (flow >= impeded_max) & (flow < diffuse_max)] = 2
    classes[valid & (flow >= diffuse_max) & (flow <= channelized_max)] = 3
    classes[valid & (flow > channelized_max)] = 4

    if not set(np.unique(classes[valid])).issubset({1, 2, 3, 4}):
        raise ValueError('Valid cells contain unexpected class codes')

    class_counts = {
        str(code): int(np.sum(classes[valid] == code)) for code in (1, 2, 3, 4)
    }

    return classes, {
        'spatial_mean': mu,
        'spatial_population_sd': sigma,
        'impeded_upper': impeded_max,
        'diffuse_upper': diffuse_max,
        'channelized_upper': channelized_max,
        'class_counts': class_counts,
    }


def process_group(index, target, ssp, radius, overwrite):
    member_list = members()
    n_members = len(member_list)
    if n_members != 180:
        raise ValueError(f'Expected 180 members; found {n_members}')

    output_dir = OUTPUT_ROOT / f'{ssp}_rad{radius}'
    output_dir.mkdir(parents=True, exist_ok=True)
    output = output_dir / \
        f'connectivity_classes_by_taxon_{ssp}_rad{radius}.npz'
    metadata_path = (
        output_dir /
        f'connectivity_classes_by_taxon_metadata_{ssp}_rad{radius}.json'
    )

    if output.is_file() and not overwrite:
        print(f'SKIP existing {output}', flush=True)
        return

    check_completeness(index, ssp, radius)

    surfaces = [
        member_surface(index, member, ssp, radius, target)
        for member in member_list
    ]
    flow_stack = np.stack(surfaces)
    common = target['footprint'] & np.all(np.isfinite(flow_stack), axis=0)
    if not common.any():
        raise ValueError(f'No common valid footprint: {ssp}, rad{radius}')
    del flow_stack

    class_stack = []
    records = []

    for position, (member, flow) in enumerate(zip(member_list, surfaces)):
        classes, stats = classify_member(flow, common)
        class_stack.append(classes)
        records.append({
            'member_index': position,
            'member_id': member_id(member),
            'pooled_member_id': pooled_member_id(member),
            'taxon': member[0],
            'sdm': member[1],
            'gcm': member[2],
            'lum': member[3],
            'curve': member[4],
            **stats,
        })
        print(
            f'{ssp} rad{radius} [{position + 1:03d}/{n_members}] '
            f'{member_id(member)}',
            flush=True,
        )

    class_stack = np.stack(class_stack).astype(np.uint8)
    expected_shape = (n_members, *target['shape'])
    if class_stack.shape != expected_shape:
        raise ValueError(
            f'Expected {expected_shape}; found {class_stack.shape}')

    ids = np.asarray([member_id(m) for m in member_list])

    temporary = output.with_name(output.stem + '.tmp.npz')
    np.savez_compressed(
        temporary,
        connectivity_class=class_stack,
        member_id=ids,
        pooled_member_id=np.asarray([pooled_member_id(m)
                                    for m in member_list]),
        taxon=np.asarray([m[0] for m in member_list]),
        sdm=np.asarray([m[1] for m in member_list]),
        gcm=np.asarray([m[2] for m in member_list]),
        lum=np.asarray([m[3] for m in member_list]),
        curve=np.asarray([m[4] for m in member_list]),
        ssp=np.asarray([ssp] * n_members),
        radius=np.asarray([radius] * n_members, dtype=np.int16),
        common_valid=common.astype(np.uint8),
        nodata_class=np.asarray(NODATA_CLASS, dtype=np.uint8),
    )
    temporary.replace(output)

    # Verify the final archive immediately.
    with np.load(output, allow_pickle=False) as archive:
        saved_classes = archive['connectivity_class']
        saved_ids = archive['member_id'].astype(str)
        saved_taxa = archive['taxon'].astype(str)

    if saved_classes.shape != expected_shape:
        raise ValueError('Saved connectivity archive shape validation failed')
    if not np.array_equal(saved_ids, ids):
        raise ValueError('Saved member order validation failed')
    for taxon in TAXA:
        if int(np.sum(saved_taxa == taxon)) != n_members // len(TAXA):
            raise ValueError(f'Unexpected member count for taxon {taxon}')

    metadata = {
        'ssp': ssp,
        'radius_pixels': radius,
        'radius_km': radius * 50,
        'member_count': n_members,
        'members_per_taxon': n_members // len(TAXA),
        'taxa': TAXA,
        'class_codes': {
            '1': 'impeded',
            '2': 'diffuse',
            '3': 'channelized',
            '4': 'intensified',
            '255': 'NoData',
        },
        'threshold_scope': 'taxon-specific member',
        'member_id_format': 'taxon__sdm__gcm__lum__curve',
        'pooled_member_id_format': 'sdm__gcm__lum__curve',
        'common_valid_footprint': 'finite in all 180 taxon-specific surfaces',
        'common_valid_cells': int(common.sum()),
        'output': str(output),
        'members': records,
    }
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding='utf-8')

    print(f'WROTE {output}', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ssp', nargs='+', choices=SSPS, default=SSPS)
    parser.add_argument('--radius', nargs='+', type=int,
                        choices=RADII, default=RADII)
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()

    target = load_template()
    index = discover()

    for ssp in args.ssp:
        for radius in args.radius:
            process_group(index, target, ssp, radius, args.overwrite)


if __name__ == '__main__':
    main()
