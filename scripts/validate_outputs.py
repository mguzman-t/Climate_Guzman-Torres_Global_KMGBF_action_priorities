#!/usr/bin/env python3
import json

import numpy as np
import rasterio
from config import (
    ANALYSIS_ROOT, SSPS, RADII, ACTION_CODES, TIE_CODE, EXPECTED_MEMBER_COUNT,
)

REQUIRED = [
    'kmgbf_modal_priority.tif',
    'kmgbf_modal_priority_agreement_percent.tif',
    'valid_realization_count.tif',
    'connectivity_impeded_agreement_percent.tif',
    'connectivity_diffuse_agreement_percent.tif',
    'connectivity_channelized_agreement_percent.tif',
    'connectivity_intensified_agreement_percent.tif',
    'connectivity_concentrated_agreement_percent.tif',
    'target2_mean_recovery_score.tif',
    'target2_mean_recovery_alignment_index.tif',
    'target2_mean_restoration_gap_index.tif',
    'metadata.json',
]


def main():
    failures = []
    reference = None
    allowed = set([0, *ACTION_CODES, TIE_CODE])

    for ssp in SSPS:
        for rad in RADII:
            folder = ANALYSIS_ROOT / f'{ssp}_rad{rad}'
            print(f'CHECK {ssp} rad{rad}')

            for name in REQUIRED:
                if not (folder / name).is_file():
                    failures.append(f'Missing {folder / name}')

            for code in ACTION_CODES:
                for name in [f'policy_selection_percent_code{code}.tif',
                             f'policy_mean_normalized_score_code{code}.tif']:
                    if not (folder / name).is_file():
                        failures.append(f'Missing {folder / name}')

            # Ensemble size from metadata.
            meta_path = folder / 'metadata.json'
            if meta_path.is_file():
                meta = json.loads(meta_path.read_text(encoding='utf-8'))
                n = int(meta.get('member_count', -1))
                print(
                    f'  members: {n} | taxa: {meta.get("connectivity_taxa")}')
                if n != EXPECTED_MEMBER_COUNT:
                    failures.append(
                        f'member_count {n} != {EXPECTED_MEMBER_COUNT} in {meta_path}')

            # Valid member count must not exceed ensemble size.
            count_path = folder / 'valid_realization_count.tif'
            if count_path.is_file():
                with rasterio.open(count_path) as src:
                    counts = src.read(1, masked=True).compressed()
                if counts.size:
                    cmin, cmax = int(counts.min()), int(counts.max())
                    print(f'  valid members per cell: min={cmin} max={cmax}')
                    if cmax > EXPECTED_MEMBER_COUNT:
                        failures.append(
                            f'Valid count {cmax} exceeds '
                            f'{EXPECTED_MEMBER_COUNT} in {count_path}')

            modal = folder / 'kmgbf_modal_priority.tif'
            if modal.is_file():
                with rasterio.open(modal) as src:
                    sig = (src.crs, src.transform, src.shape)
                    if reference is None:
                        reference = sig
                    elif sig != reference:
                        failures.append(f'Grid mismatch: {modal}')
                    vals = src.read(1, masked=True).compressed()
                bad = set(np.unique(vals).astype(int)) - allowed
                if bad:
                    failures.append(
                        f'Unexpected codes {sorted(bad)} in {modal}')
                u, n = np.unique(vals, return_counts=True)
                print('  modal:', ' '.join(
                    f'{int(c)}={100 * k / len(vals):.2f}%' for c, k in zip(u, n)))

    if failures:
        print('\nFAILED')
        for x in failures:
            print(' -', x)
        raise SystemExit(1)
    print('\nALL CHECKS PASSED')


if __name__ == '__main__':
    main()
