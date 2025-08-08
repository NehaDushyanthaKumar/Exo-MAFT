"""
FeatureTracker module for plotting and exporting molecular absorption features
on a combined transmission spectrum using a HITRAN-style .par linelist.

This module provides:
  - Load and parse transmission spectra
  - Load and process HITRAN linelists
  - Plot feature lines overlaid on a spectrum
  - Save output plots and line lists

Output:
  - plots/feature_tracking.png
  - output/feature_lines.txt
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class FeatureTracker:
    """
    Class to track and visualize molecular spectral features in transmission spectra.

    Attributes
    ----------
    combined_spectrum_path : str
        Path to the combined spectrum data file.
    par_file_path : str
        Path to the HITRAN-style .par linelist file.
    plot_file : str
        Path to save the feature plot image.
    out_lines_file : str
        Path to save the molecular feature lines.
    molecules : dict
        Dictionary of molecules with HITRAN IDs.
    top_n : int
        Number of strongest lines to select per molecule.
    bin_width : float
        Bin width for clustering lines in microns.
    alpha_lines : float
        Transparency of the vertical feature lines on the plot.
    """

    def __init__(self):
        self.combined_spectrum_path = 'output/combined_spectrum.txt'
        self.par_file_path = 'data/6894c8ca.par'
        self.plot_file = 'plots/feature_tracking.png'
        self.out_lines_file = 'output/feature_lines.txt'

        self.molecules = {
            'H2O':  1,
            'CO2':  2,
            'NH3':  3,
            'CO':   5,
            'CH4':  6,
            'HCN':  8,
        }

        self.top_n = 30
        self.bin_width = 0.02
        self.alpha_lines = 0.3

        os.makedirs('plots', exist_ok=True)
        os.makedirs('output', exist_ok=True)

    def load_combined_spectrum(self):
        """Load the transmission spectrum from a tab-separated text file."""
        return pd.read_csv(
            self.combined_spectrum_path,
            sep='\t',
            skiprows=1,
            names=['wavelength_um','depth','depth_err'],
            comment='#',
            dtype={'wavelength_um':float,'depth':float,'depth_err':float}
        )

    def load_linelist_par(self):
        """
        Load a HITRAN .par linelist file using fixed-width formatting.

        Returns
        -------
        pandas.DataFrame
            DataFrame containing molecular lines and converted wavelengths.
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
            self.par_file_path,
            colspecs=colspecs,
            names=names,
            comment='#',
            skip_blank_lines=True
        )

        for c in ('molec_id','local_iso_id'):
            bad_rows = df[~df[c].astype(str).str.fullmatch(r'\d+')]
            if not bad_rows.empty:
                print(f"\nColumn {c} has non-integer values:")
                print(bad_rows[[c]])
            else:
                df[c] = df[c].astype(float)

        df['wavelength_um'] = 1e4 / df['nu']
        return df

    def plot_and_save(self):
        """
        Create the feature line plot, save the figure and output the feature lines.
        """
        spec = self.load_combined_spectrum()
        lines = self.load_linelist_par()

        wave = spec['wavelength_um'].values
        depth = spec['depth'].values

        fig, ax = plt.subplots(figsize=(10, 4))
        ax.plot(wave, depth, color='k', lw=1, label='Combined spectrum')
        ax.fill_between(wave, depth - spec['depth_err'], depth + spec['depth_err'],
                        color='gray', alpha=0.2)

        ymin, ymax = ax.get_ylim()
        out_rows = []
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

        for idx, (name, mol_id) in enumerate(self.molecules.items()):
            dfm = lines[lines['molec_id'] == mol_id]
            if dfm.empty:
                continue

            strongest = dfm.nlargest(self.top_n, 'sw')['wavelength_um'].values
            bins = np.arange(wave.min(), wave.max() + self.bin_width, self.bin_width)
            inds = np.digitize(strongest, bins)
            centers = [strongest[inds == b].mean() for b in np.unique(inds)]

            for lam in centers:
                out_rows.append((name, lam))

            for lam in centers:
                ax.vlines(lam, ymin, ymax, color=colors[idx % len(colors)],
                          linewidth=1, alpha=self.alpha_lines)

            ax.plot([], [], color=colors[idx % len(colors)], lw=3,
                    alpha=self.alpha_lines, label=name)

        df_out = pd.DataFrame(out_rows, columns=['molecule', 'wavelength_um'])
        df_out.to_csv(self.out_lines_file, sep='\t', index=False, float_format='%.6f')

        ax.set_xlim(wave.min(), wave.max())
        ax.set_ylim(ymin, ymax)
        ax.set_xlabel('Wavelength (\u03bcm)')
        ax.set_ylabel('Transit Depth (Rp\u00b2/Rs\u00b2)')
        ax.set_title('Combined Spectrum with Molecular Feature Lines')
        ax.legend(loc='upper right', frameon=False, fontsize='small', ncol=3)
        ax.grid(alpha=0.2)

        plt.tight_layout()
        plt.savefig(self.plot_file, dpi=200)
        plt.show()


if __name__ == "__main__":
    tracker = FeatureTracker()
    tracker.plot_and_save()
