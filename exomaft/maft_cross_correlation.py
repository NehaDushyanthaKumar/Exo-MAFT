#!/usr/bin/env python3
"""
cross_correlation.py

This module defines the CrossCorrelator class to compute cross-correlation between a combined transmission
spectrum and molecular line templates. It generates cross-correlation function plots and outputs peak values.

The class structure improves modularity and reusability.

Outputs:
  • plots/ccf_all.png
  • output/ccf_peaks.txt
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import correlate


class CrossCorrelator:
    """
    Class to perform cross-correlation between a combined transmission spectrum
    and molecular line templates.

    Attributes
    ----------
    combined_spectrum : str
        File path to the combined spectrum data.
    par_file : str
        File path to the molecular line parameter file.
    molecules : dict
        Dictionary of molecule names and their corresponding line IDs.
    top_n : int
        Number of top strongest lines to use per molecule.
    bin_width : float
        Width of the wavelength bins used when clustering strongest lines.
    plot_file : str
        File path for saving the cross-correlation function plot.
    out_peaks : str
        File path for saving cross-correlation peak values.
    """

    def __init__(self,
                 combined_spectrum='exomaft/output/combined_spectrum.txt',
                 par_file='exomaft/data/6894c8ca.par',
                 plot_file='plots/ccf_all.png',
                 out_peaks='output/ccf_peaks.txt',
                 molecules=None,
                 top_n=30,
                 bin_width=0.02):
        """
        Initialize the CrossCorrelator instance.

        Parameters
        ----------
        combined_spectrum : str, optional
            Path to the combined spectrum file (default is 'exomaft/output/combined_spectrum.txt').
        par_file : str, optional
            Path to the molecular line list parameter file (default is 'exomaft/data/6894c8ca.par').
        plot_file : str, optional
            Path to save the cross-correlation function plot (default is 'plots/ccf_all.png').
        out_peaks : str, optional
            Path to save the cross-correlation peaks text file (default is 'output/ccf_peaks.txt').
        molecules : dict, optional
            Dictionary mapping molecule names to molecular ID (default uses a predefined set).
        top_n : int, optional
            Number of top strongest lines to consider per molecule (default is 30).
        bin_width : float, optional
            Wavelength bin width used to cluster lines (default is 0.02 µm).
        """
        self.combined_spectrum = combined_spectrum
        self.par_file = par_file
        self.plot_file = plot_file
        self.out_peaks = out_peaks
        self.molecules = molecules or {
            'H2O': 1,
            'CO2': 2,
            'NH3': 3,
            'CO':  5,
            'CH4': 6,
            'HCN': 8,
        }
        self.top_n = top_n
        self.bin_width = bin_width

        # Ensure output directories exist
        os.makedirs(os.path.dirname(self.plot_file), exist_ok=True)
        os.makedirs(os.path.dirname(self.out_peaks), exist_ok=True)

    def load_combined_spectrum(self):
        """
        Load the combined transmission spectrum data.

        Returns
        -------
        tuple of (np.ndarray, np.ndarray)
            wavelength array and depth array
        """
        spec = pd.read_csv(
            self.combined_spectrum,
            sep='\t',
            skiprows=1,
            names=['wavelength_um', 'depth', 'depth_err'],
            comment='#',
            dtype={'wavelength_um': float, 'depth': float, 'depth_err': float}
        )
        return spec['wavelength_um'].values, spec['depth'].values

    def load_linelist_par(self):
        """
        Load the molecular line list parameter file.

        Returns
        -------
        pandas.DataFrame
            DataFrame with molecular line parameters including wavelength_um column.
        """
        colspecs = [(0, 2), (2, 3), (3, 15), (15, 25), (25, 35), (35, 40),
                    (40, 45), (45, 55), (55, 59), (59, 67)]
        names = ['molec_id', 'local_iso_id', 'nu', 'sw', 'a', 'gamma_air',
                 'gamma_self', 'elower', 'n_air', 'delta_air']
        df = pd.read_fwf(self.par_file, colspecs=colspecs, names=names,
                         comment='#', skip_blank_lines=True)
        # Correct data types
        df['molec_id'] = df['molec_id'].astype(int)
        df['local_iso_id'] = df['local_iso_id'].astype(int)
        for col in ['nu', 'sw', 'a', 'gamma_air', 'gamma_self', 'elower', 'n_air', 'delta_air']:
            df[col] = df[col].astype(float)
        # Convert wavenumber to wavelength in microns
        df['wavelength_um'] = 1e4 / df['nu']
        return df

    @staticmethod
    def zscore(arr):
        """
        Compute the z-score normalization of the input array.

        Parameters
        ----------
        arr : np.ndarray
            Input array.

        Returns
        -------
        np.ndarray
            Z-scored array.
        """
        std = np.std(arr)
        return (arr - np.mean(arr)) / (std if std > 0 else 1.0)

    def run_cross_correlation(self):
        """
        Perform cross-correlation between the combined spectrum and molecular templates.
        Creates plot and saves peak correlation values.
        """
        wave, depth = self.load_combined_spectrum()
        lines = self.load_linelist_par()

        ccf_peaks = []
        plt.figure(figsize=(10, 4))

        for name, mol_id in self.molecules.items():
            dfm = lines[lines['molec_id'] == mol_id]
            if dfm.empty:
                continue

            # Select strongest lines for clustering
            strongest = dfm.nlargest(self.top_n, 'sw')['wavelength_um'].values
            bins = np.arange(wave.min(), wave.max() + self.bin_width, self.bin_width)
            inds = np.digitize(strongest, bins)
            centers = np.array([strongest[inds == b].mean() for b in np.unique(inds)])

            # Create template vector with spikes at cluster centers
            tpl = np.zeros_like(depth)
            for lam in centers:
                i = np.argmin(np.abs(wave - lam))
                tpl[i] = 1.0

            # Compute z-scores
            d0 = self.zscore(depth)
            tpl0 = self.zscore(tpl)

            # Full-mode correlation
            ccf = correlate(d0, tpl0, mode='same')
            peak = np.max(ccf)
            ccf_peaks.append((name, peak))

            # Plot CCF trace for this molecule
            plt.plot(wave, ccf, label=f"{name} (peak {peak:.1f})", linewidth=1)

        # Save peak values
        df_peaks = pd.DataFrame(ccf_peaks, columns=['molecule', 'ccf_peak'])
        df_peaks.to_csv(self.out_peaks, sep='\t', index=False, float_format='%.3f')

        # Finalize plot
        plt.xlabel('Wavelength (µm)')
        plt.ylabel('CCF (z-score correlation)')
        plt.title('Cross-Correlation Functions (All Molecules)')
        plt.legend(loc='upper right', fontsize='small', ncol=2, frameon=False)
        plt.grid(alpha=0.2)
        plt.tight_layout()
        plt.savefig(self.plot_file, dpi=200)
        plt.show()


if __name__ == '__main__':
    correlator = CrossCorrelator()
    correlator.run_cross_correlation()
