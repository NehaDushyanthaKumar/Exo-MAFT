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

from pathlib import Path
HERE = Path(__file__).resolve().parent        # .../Exo-MAFT/exomaft
ROOT = HERE.parent                            # .../Exo-MAFT

DATA_DIR = HERE / 'data'
PAR_FILE = DATA_DIR / '6894c8ca.par'          # <-- no extra 'exomaft'

PLOT_DIR = ROOT / 'plots'                     # existing folders at repo root
OUT_DIR  = ROOT / 'output'

COMBINED_SPECTRUM = OUT_DIR / 'combined_spectrum.txt'  # change if yours is elsewhere
PLOT_FILE = PLOT_DIR / 'feature_tracking.png'
OUT_LINES = OUT_DIR  / 'feature_lines.txt'

# sanity checks (fail fast if wrong)
for p in [PAR_FILE, COMBINED_SPECTRUM]:
    if not p.exists():
        raise FileNotFoundError(f"Expected file not found: {p}")

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


print("Reading PAR_FILE:", PAR_FILE)
with open(PAR_FILE, 'rb') as f:
    first = f.readline()
print("First 30 bytes:", first[:30])

# =============================================================================
# 3) LOAD HITRAN‐STYLE LINELIST (.par fixed-width)
# =============================================================================
def load_linelist_par(path):
    colspecs = [
        (0,2),(2,3),(3,15),(15,25),(25,35),
        (35,40),(40,45),(45,55),(55,59),(59,67)
    ]
    names = [
        'molec_id','local_iso_id','nu','sw','a',
        'gamma_air','gamma_self','elower','n_air','delta_air'
    ]

    # Read as strings to avoid crashing on odd lines; skip blank/comment lines.
    df = pd.read_fwf(
        path,
        colspecs=colspecs,
        names=names,
        dtype=str,
        comment='#',
        skip_blank_lines=True,
        on_bad_lines='skip'  # pandas>=1.3
    )

    # Keep only rows where IDs are clean integers
    is_int = lambda s: s.str.strip().str.fullmatch(r'\d+').fillna(False)
    df = df[is_int(df['molec_id']) & is_int(df['local_iso_id'])].copy()

    # Cast IDs to int, others to numeric
    df[['molec_id','local_iso_id']] = df[['molec_id','local_iso_id']].astype(int)
    for c in ['nu','sw','a','gamma_air','gamma_self','elower','n_air','delta_air']:
        df[c] = pd.to_numeric(df[c], errors='coerce')

    # Require at least nu & sw
    df = df.dropna(subset=['nu','sw'])
    df['wavelength_um'] = 1e4 / df['nu']
    return df

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
