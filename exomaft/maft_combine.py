#!/usr/bin/env python3
"""
process_and_plot_spectra.py

This script reads multiple transmission spectrum files (NIRISS, Combined Model,
Archival HST/Spitzer, and NIRSpec PRISM), standardizes and merges them into a
single combined spectrum (removing overlapping wavelengths), and generates
plots (with and without error bars) saved to disk.

Outputs:
  - output/combined_spectrum.txt          : merged spectrum file
  - plots/overlay_spectra_with_errors.png : overlay plot with 1σ error bars
  - plots/overlay_spectra_no_errors.png   : overlay plot without error bars
  - plots/combined_spectrum.png           : combined spectrum plot with error handling

Usage:
  Update the input file paths as needed and run:
    python process_and_plot_spectra.py
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# 1. CONFIGURE INPUT & OUTPUT PATHS
# =============================================================================
# Input spectrum files (update paths as necessary)
file_niriss = '/Users/nehadushyanthakumar/Desktop/FIREFLy_transmission_spectrum_R_300_order_1.txt'
file_comb   = '/Users/nehadushyanthakumar/Desktop/Exo-MAFT/data/ZENODO/5_MODEL_SPECTRUM/BestFitModel.txt'
file_arch   = '/Users/nehadushyanthakumar/Desktop/Exo-MAFT/data/ZENODO/4_TRANSMISSION_SPECTRA/Archival_Data.txt'
file_prism  = '/Users/nehadushyanthakumar/Desktop/FIREFLy_transit_spectrum.txt'

# Output directories
output_dir = 'output'
plots_dir  = 'plots'
os.makedirs(output_dir, exist_ok=True)
eos.makedirs(plots_dir, exist_ok=True)

# =============================================================================
# 2. LOAD AND PREPARE INDIVIDUAL SPECTRA
# =============================================================================
# 2.1 NIRISS SOSS (stored in ppm; convert to fraction)
df_niriss = pd.read_csv(
    file_niriss,
    delim_whitespace=True,
    comment='#',
    header=None,
    names=['wavelength','wavelength_err','depth_ppm','depth_err_ppm']
)
df_niriss['depth']     = df_niriss['depth_ppm'] / 1e6
(df_niriss['depth_err']) = df_niriss['depth_err_ppm'] / 1e6

# 2.2 Combined Multi-Instrument model (no intrinsic errors)
df_comb = pd.read_csv(
    file_comb,
    delim_whitespace=True,
    comment='#',
    header=None,
    names=['wavelength','depth']
)
df_comb['wavelength_err'] = 0.0
df_comb['depth_err']      = 0.0

# 2.3 Archival HST & Spitzer (WASP-39b)
df_arch = pd.read_csv(
    file_arch,
    delim_whitespace=True,
    comment='#',
    header=None,
    names=['wavelength','wavelength_err','depth','depth_err','ntransits']
)[['wavelength','wavelength_err','depth','depth_err']]

# 2.4 NIRSpec PRISM (no wavelength error provided)
df_prism = pd.read_csv(
    file_prism,
    delim_whitespace=True,
    comment='#',
    header=None,
    names=['wavelength','depth','depth_err']
)
df_prism['wavelength_err'] = 0.0

# =============================================================================
# 3. ALIGN BASELINES FOR VISUAL COMPARISON
# =============================================================================
# Shift each spectrum so that their median depth matches the NIRISS median
base_med = df_niriss['depth'].median()
for df in (df_comb, df_arch, df_prism):
    df['depth'] += (base_med - df['depth'].median())

# =============================================================================
# 4. MERGE SPECTRA INTO ONE COMBINED TABLE
# =============================================================================
# Standardize blocks
blocks = []
for df in (df_niriss, df_arch, df_prism, df_comb):
    tmp = df[['wavelength','depth']].copy()
    # preserve error if present
    tmp['depth_err'] = df.get('depth_err', np.nan)
    blocks.append(tmp)
# Concatenate, sort, remove duplicate wavelengths (keep first)
combined = pd.concat(blocks, ignore_index=True)
combined = combined.sort_values('wavelength')
combined = combined.drop_duplicates(subset='wavelength', keep='first')
# Write combined spectrum to file
out_file = os.path.join(output_dir, 'combined_spectrum.txt')
combined.to_csv(
    out_file,
    sep='\t',
    index=False,
    header=['wavelength','depth','depth_err'],
    float_format='%.6e',
    na_rep=''  # blank for missing errors
)
print(f"✅ Combined spectrum written to {out_file}")

# =============================================================================
# 5. PLOT OVERLAY WITH ERROR BARS
# =============================================================================
plt.figure(figsize=(8,5))
# NIRISS
plt.errorbar(
    df_niriss['wavelength'], df_niriss['depth'],
    xerr=df_niriss['wavelength_err'], yerr=df_niriss['depth_err'],
    fmt='o-', label='NIRISS SOSS', color='C0', capsize=2
)
# Combined
plt.errorbar(
    df_comb['wavelength'], df_comb['depth'],
    xerr=df_comb['wavelength_err'], yerr=df_comb['depth_err'],
    fmt='s--', label='Combined Model', color='C1', capsize=2
)
# Archival
plt.errorbar(
    df_arch['wavelength'], df_arch['depth'],
    xerr=df_arch['wavelength_err'], yerr=df_arch['depth_err'],
    fmt='^:', label='Archival HST/Spitzer', color='C2', capsize=2
)
# PRISM
plt.errorbar(
    df_prism['wavelength'], df_prism['depth'],
    xerr=df_prism['wavelength_err'], yerr=df_prism['depth_err'],
    fmt='D-.', label='NIRSpec PRISM', color='C3', capsize=2
)
plt.xlabel('Wavelength (μm)')
plt.ylabel('Transit Depth (Rp²/Rs²)')
plt.title('Overlay of Spectra with 1σ Error Bars')
plt.legend(loc='best', fontsize='small')
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(plots_dir, 'overlay_spectra_with_errors.png'), dpi=200)
plt.close()
print(f"✅ Saved plot with errors to {os.path.join(plots_dir, 'overlay_spectra_with_errors.png')}")

# =============================================================================
# 6. PLOT OVERLAY WITHOUT ERROR BARS
# =============================================================================
plt.figure(figsize=(8,5))
plt.plot(df_niriss['wavelength'], df_niriss['depth'], '-o', label='NIRISS SOSS', color='C0')
plt.plot(df_comb['wavelength'], df_comb['depth'], '--s', label='Combined Model', color='C1')
plt.plot(df_arch['wavelength'], df_arch['depth'], '-.^', label='Archival HST/Spitzer', color='C2')
plt.plot(df_prism['wavelength'], df_prism['depth'], '-.D', label='NIRSpec PRISM', color='C3')
plt.xlabel('Wavelength (μm)')
plt.ylabel('Transit Depth (Rp²/Rs²)')
plt.title('Overlay of Spectra (No Error Bars)')
plt.legend(loc='best', fontsize='small')
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(plots_dir, 'overlay_spectra_no_errors.png'), dpi=200)
plt.close()
print(f"✅ Saved plot without errors to {os.path.join(plots_dir, 'overlay_spectra_no_errors.png')}")

# =============================================================================
# 7. PLOT COMBINED SPECTRUM
# =============================================================================
plt.figure(figsize=(8,5))
# plot points with errors
with_err = combined[combined['depth_err'].notna()]
no_err   = combined[combined['depth_err'].isna()]
if not with_err.empty:
    plt.errorbar(
        with_err['wavelength'], with_err['depth'], yerr=with_err['depth_err'],
        fmt='o', label='With error bars', color='black', capsize=2
    )
# plot points without errors
if not no_err.empty:
    plt.plot(
        no_err['wavelength'], no_err['depth'], 's',
        label='No error bar', color='blue'
    )
plt.xlabel('Wavelength (μm)')
plt.ylabel('Transit Depth (Rp²/Rs²)')
plt.title('Combined Transmission Spectrum')
plt.legend(loc='best', fontsize='small')
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(plots_dir, 'combined_spectrum.png'), dpi=200)
plt.close()
print(f"✅ Saved combined spectrum plot to {os.path.join(plots_dir, 'combined_spectrum.png')}")
