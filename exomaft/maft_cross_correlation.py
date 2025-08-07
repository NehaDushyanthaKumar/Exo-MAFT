#!/usr/bin/env python3
"""
cross_correlation.py

Compute cross-correlation between the combined transmission spectrum
and molecular line templates, then save:
  • plots/ccf_all.png
  • output/ccf_peaks.txt
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import correlate

# =============================================================================
# CONFIGURATION
# =============================================================================
COMBINED_SPECTRUM = 'data/combined_spectrum.txt'
PAR_FILE          = 'data/6894c8ca.par'
PLOT_FILE         = 'plots/ccf_all.png'
OUT_PEAKS         = 'output/ccf_peaks.txt'

MOLECULES = {
    'H2O': 1,
    'CO2': 2,
    'NH3': 3,
    'CO':  5,
    'CH4': 6,
    'HCN': 8,
}

TOP_N     = 30
BIN_WIDTH = 0.02

# =============================================================================
# 1) MAKE OUTPUT FOLDERS
# =============================================================================
os.makedirs('plots',  exist_ok=True)
os.makedirs('output', exist_ok=True)

# =============================================================================
# 2) LOAD SPECTRUM
# =============================================================================
spec = pd.read_csv(
    COMBINED_SPECTRUM,
    sep='\t',
    skiprows=1,
    names=['wavelength_um','depth','depth_err'],
    comment='#',
    dtype={'wavelength_um':float,'depth':float,'depth_err':float}
)
wave  = spec['wavelength_um'].values
depth = spec['depth'].values

# =============================================================================
# 3) LOAD LINELIST (same as above)
# =============================================================================
def load_linelist_par(path):
    colspecs = [(0,2),(2,3),(3,15),(15,25),(25,35),(35,40),(40,45),(45,55),(55,59),(59,67)]
    names    = ['molec_id','local_iso_id','nu','sw','a','gamma_air','gamma_self','elower','n_air','delta_air']
    df = pd.read_fwf(path, colspecs=colspecs, names=names, comment='#', skip_blank_lines=True)
    for c in names:
        if c in ('molec_id','local_iso_id'):
            df[c] = df[c].astype(int)
        else:
            df[c] = df[c].astype(float)
    df['wavelength_um'] = 1e4 / df['nu']
    return df

lines = load_linelist_par(PAR_FILE)

# =============================================================================
# 4) Z‐scoring UTILITY
# =============================================================================
def zscore(arr):
    """Return (arr − mean)/std, guarding against zero std."""
    std = np.std(arr)
    return (arr - np.mean(arr)) / (std if std>0 else 1.0)

# =============================================================================
# 5) CROSS‐CORRELATE & RECORD PEAKS
# =============================================================================
ccf_peaks = []

plt.figure(figsize=(10,4))

for idx, (name, mol_id) in enumerate(MOLECULES.items()):
    dfm = lines[lines['molec_id']==mol_id]
    if dfm.empty:
        continue

    # strongest lines → cluster
    strongest = dfm.nlargest(TOP_N, 'sw')['wavelength_um'].values
    bins      = np.arange(wave.min(), wave.max()+BIN_WIDTH, BIN_WIDTH)
    inds      = np.digitize(strongest, bins)
    centers   = np.array([ strongest[inds==b].mean() for b in np.unique(inds) ])

    # template: zeros everywhere except spikes at 'centers'
    tpl = np.zeros_like(depth)
    for lam in centers:
        i = np.argmin(np.abs(wave - lam))
        tpl[i] = 1.0

    # compute z‐scores
    d0   = zscore(depth)
    tpl0 = zscore(tpl)

    # full‐mode correlation
    ccf = correlate(d0, tpl0, mode='same')
    peak = np.max(ccf)
    ccf_peaks.append((name, peak))

    # plot each CCF trace
    plt.plot(wave, ccf, label=f"{name} (peak {peak:.1f})", linewidth=1)

# =============================================================================
# 6) SAVE PEAKS TO TEXT
# =============================================================================
df_peaks = pd.DataFrame(ccf_peaks, columns=['molecule','ccf_peak'])
df_peaks.to_csv(
    OUT_PEAKS,
    sep='\t',
    index=False,
    float_format='%.3f'
)

# =============================================================================
# 7) FINALIZE & SAVE FIGURE
# =============================================================================
plt.xlabel('Wavelength (µm)')
plt.ylabel('CCF (z-score correlation)')
plt.title('Cross-Correlation Functions (All Molecules)')
plt.legend(loc='upper right', fontsize='small', ncol=2, frameon=False)
plt.grid(alpha=0.2)
plt.tight_layout()
plt.savefig(PLOT_FILE, dpi=200)
plt.show()
