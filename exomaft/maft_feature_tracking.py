#!/usr/bin/env python3
"""
feature_tracking.py

Load a combined transmission spectrum and a HITRAN‐style .par linelist,
then overlay semi‐transparent molecular feature lines and save:
  • plots/feature_tracking.png
  • output/feature_lines.txt
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# =============================================================================
# CONFIGURATION
# =============================================================================
COMBINED_SPECTRUM = 'data/combined_spectrum.txt'  # input spectrum
PAR_FILE          = 'data/6894c8ca.par'          # input linelist
PLOT_FILE         = 'plots/feature_tracking.png'
OUT_LINES         = 'output/feature_lines.txt'

# Pick molecules & HITRAN IDs to track
MOLECULES = {
    'H2O':  1,
    'CO2':  2,
    'NH3':  3,
    'CO':   5,
    'CH4':  6,
    'HCN':  8,
}

TOP_N       = 30    # strongest lines per molecule
BIN_WIDTH   = 0.02  # μm -- cluster width
ALPHA_LINES = 0.3   # line transparency

# =============================================================================
# 1) MAKE OUTPUT DIRS
# =============================================================================
os.makedirs('plots',  exist_ok=True)
os.makedirs('output', exist_ok=True)

# =============================================================================
# 2) LOAD COMBINED SPECTRUM
# =============================================================================
#    skip the header row and force float conversion
spec = pd.read_csv(
    COMBINED_SPECTRUM,
    sep='\t',
    skiprows=1,
    names=['wavelength_um','depth','depth_err'],
    comment='#',
    dtype={'wavelength_um':float,'depth':float,'depth_err':float}
)
wave  = spec['wavelength_um'].values  # array of wavelengths
depth = spec['depth'].values           # array of depths

# =============================================================================
# 3) LOAD HITRAN‐STYLE LINELIST (.par fixed-width)
# =============================================================================
def load_linelist_par(path):
    """
    Read a HITRAN .par file via pandas.read_fwf, assign column names,
    cast types, and compute wavelength in microns.
    """
    colspecs = [
        (0,2),(2,3),(3,15),(15,25),(25,35),
        (35,40),(40,45),(45,55),(55,59),(59,67)
    ]
    names = [
        'molec_id','local_iso_id','nu','sw','a',
        'gamma_air','gamma_self','elower','n_air','delta_air'
    ]
    df = pd.read_fwf(
        path,
        colspecs=colspecs,
        names=names,
        comment='#',
        skip_blank_lines=True
    )
    # cast columns to numeric
    for c in names:
        if c in ('molec_id','local_iso_id'):
            df[c] = df[c].astype(int)
        else:
            df[c] = df[c].astype(float)
    # convert wavenumber (cm⁻¹) → wavelength (µm)
    df['wavelength_um'] = 1e4 / df['nu']
    return df

lines = load_linelist_par(PAR_FILE)

# =============================================================================
# 4) OPEN FIGURE & PLOT SPECTRUM
# =============================================================================
fig, ax = plt.subplots(figsize=(10,4))

# Draw the combined spectrum (black line) with error region
ax.plot(wave, depth, color='k', lw=1, label='Combined spectrum')
ax.fill_between(
    wave,
    depth - spec['depth_err'],
    depth + spec['depth_err'],
    color='gray', alpha=0.2
)

# get vertical extent for feature lines
ymin, ymax = ax.get_ylim()

# =============================================================================
# 5) CLUSTER & OVERLAY MOLECULAR LINES
# =============================================================================
# prepare an output table
out_rows = []

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
for idx, (name, mol_id) in enumerate(MOLECULES.items()):
    # select lines for this molecule
    dfm = lines[lines['molec_id'] == mol_id]
    if dfm.empty:
        continue

    # pick top‐N by line strength 'sw'
    strongest = dfm.nlargest(TOP_N, 'sw')['wavelength_um'].values

    # cluster into BIN_WIDTH bins
    bins   = np.arange(wave.min(), wave.max() + BIN_WIDTH, BIN_WIDTH)
    inds   = np.digitize(strongest, bins)
    centers = [ strongest[inds == b].mean() for b in np.unique(inds) ]

    # record to output table
    for lam in centers:
        out_rows.append((name, lam))

    # draw semi‐transparent full‐height lines
    for lam in centers:
        ax.vlines(
            lam, ymin, ymax,
            color=colors[idx % len(colors)],
            linewidth=1,
            alpha=ALPHA_LINES
        )

    # add dummy for the legend
    ax.plot(
        [], [], color=colors[idx % len(colors)],
        lw=3, alpha=ALPHA_LINES, label=name
    )

# =============================================================================
# 6) SAVE FEATURE LIST TO .TXT
# =============================================================================
df_out = pd.DataFrame(out_rows, columns=['molecule','wavelength_um'])
df_out.to_csv(
    OUT_LINES,
    sep='\t',
    index=False,
    float_format='%.6f'
)

# =============================================================================
# 7) FINALIZE & SAVE FIGURE
# =============================================================================
ax.set_xlim(wave.min(), wave.max())
ax.set_ylim(ymin, ymax)
ax.set_xlabel('Wavelength (µm)')
ax.set_ylabel('Transit Depth (Rp²/Rs²)')
ax.set_title('Combined Spectrum with Molecular Feature Lines')
ax.legend(loc='upper right', frameon=False, fontsize='small', ncol=3)
ax.grid(alpha=0.2)
plt.tight_layout()
plt.savefig(PLOT_FILE, dpi=200)
plt.show()
