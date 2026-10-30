#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import rasterio
from config import ANALYSIS_ROOT, SSPS, RADII, ACTION_CODES, TIE_CODE

REQUIRED = [
    'kmgbf_modal_priority.tif',
    'kmgbf_modal_priority_agreement_percent.tif',
    'valid_realization_count.tif',
    'connectivity_impeded_agreement_percent.tif',
    'connectivity_diffuse_agreement_percent.tif',
    'connectivity_concentrated_agreement_percent.tif',
    'target2_mean_recovery_score.tif',
    'target2_mean_recovery_alignment_index.tif',
    'target2_mean_restoration_gap_index.tif',
]

def main():
    failures=[]; reference=None
    allowed=set([0,*ACTION_CODES,TIE_CODE])
    for ssp in SSPS:
        for rad in RADII:
            folder=ANALYSIS_ROOT/f'{ssp}_rad{rad}'
            print(f'CHECK {ssp} rad{rad}')
            for name in REQUIRED:
                path=folder/name
                if not path.is_file(): failures.append(f'Missing {path}')
            for code in ACTION_CODES:
                for name in [f'policy_selection_percent_code{code}.tif',f'policy_mean_normalized_score_code{code}.tif']:
                    if not (folder/name).is_file(): failures.append(f'Missing {folder/name}')
            modal=folder/'kmgbf_modal_priority.tif'
            if modal.is_file():
                with rasterio.open(modal) as src:
                    sig=(src.crs,src.transform,src.shape)
                    if reference is None: reference=sig
                    elif sig!=reference: failures.append(f'Grid mismatch: {modal}')
                    vals=src.read(1,masked=True).compressed()
                    bad=set(np.unique(vals).astype(int))-allowed
                    if bad: failures.append(f'Unexpected codes {sorted(bad)} in {modal}')
                    u,n=np.unique(vals,return_counts=True)
                    print('  modal:', ' '.join(f'{int(c)}={100*k/len(vals):.2f}%' for c,k in zip(u,n)))
    if failures:
        print('\nFAILED')
        for x in failures: print(' -',x)
        raise SystemExit(1)
    print('\nALL CHECKS PASSED')
if __name__=='__main__': main()
