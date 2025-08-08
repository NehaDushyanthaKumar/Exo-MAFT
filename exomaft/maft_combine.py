import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

class TransmissionSpectrumProcessor:
    """
    A class to load, process, align, merge, and plot exoplanet transmission spectra 
    from different instruments including NIRISS, HST/Spitzer, NIRSpec PRISM, and a combined model.

    Parameters
    ----------
    file_paths : dict
        Dictionary containing paths to the spectrum files with keys:
        'niriss', 'combined', 'archival', and 'prism'.
    output_dir : str, optional
        Directory where the output combined spectrum will be saved. Default is 'output'.
    plots_dir : str, optional
        Directory where the output plots will be saved. Default is 'plots'.
    """
    def __init__(self, file_paths, output_dir='output', plots_dir='plots'):
        self.file_niriss = file_paths['niriss']
        self.file_comb = file_paths['combined']
        self.file_arch = file_paths['archival']
        self.file_prism = file_paths['prism']

        self.output_dir = output_dir
        self.plots_dir = plots_dir

        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(self.plots_dir, exist_ok=True)

        # DataFrames will be loaded later
        self.df_niriss = None
        self.df_comb = None
        self.df_arch = None
        self.df_prism = None
        self.combined = None

    def load_spectra(self):
        """
        Load the spectrum data files for all instruments and convert units where necessary.
        """
        # Load NIRISS spectrum and convert ppm to fractional depth
        self.df_niriss = pd.read_csv(
            self.file_niriss,
            delim_whitespace=True,
            comment='#',
            header=None,
            names=['wavelength','wavelength_err','depth_ppm','depth_err_ppm']
        )
        self.df_niriss['depth'] = self.df_niriss['depth_ppm'] / 1e6
        self.df_niriss['depth_err'] = self.df_niriss['depth_err_ppm'] / 1e6

        # Load Combined Model spectrum without errors
        self.df_comb = pd.read_csv(
            self.file_comb,
            delim_whitespace=True,
            comment='#',
            header=None,
            names=['wavelength','depth']
        )
        self.df_comb['wavelength_err'] = 0.0
        self.df_comb['depth_err'] = 0.0

        # Load archival HST and Spitzer data (include only necessary columns)
        self.df_arch = pd.read_csv(
            self.file_arch,
            delim_whitespace=True,
            comment='#',
            header=None,
            names=['wavelength','wavelength_err','depth','depth_err','ntransits']
        )[['wavelength','wavelength_err','depth','depth_err']]

        # Load NIRSpec PRISM data
        self.df_prism = pd.read_csv(
            self.file_prism,
            delim_whitespace=True,
            comment='#',
            header=None,
            names=['wavelength','depth','depth_err']
        )
        self.df_prism['wavelength_err'] = 0.0

    def align_baselines(self):
        """
        Align the baseline (median) of all datasets to match the NIRISS spectrum.
        This is useful for visually consistent overlay plots.
        """
        base_med = self.df_niriss['depth'].median()
        for df in (self.df_comb, self.df_arch, self.df_prism):
            df['depth'] += (base_med - df['depth'].median())

    def merge_spectra(self):
        """
        Merge all spectra into a single DataFrame and export it as a tab-separated file.
        """
        blocks = []
        for df in (self.df_niriss, self.df_arch, self.df_prism, self.df_comb):
            tmp = df[['wavelength','depth']].copy()
            tmp['depth_err'] = df.get('depth_err', np.nan)  # Fill missing error as NaN
            blocks.append(tmp)

        # Concatenate all data blocks and remove duplicate wavelength entries
        self.combined = pd.concat(blocks, ignore_index=True)
        self.combined = self.combined.sort_values('wavelength')
        self.combined = self.combined.drop_duplicates(subset='wavelength', keep='first')

        out_file = os.path.join(self.output_dir, 'combined_spectrum.txt')
        self.combined.to_csv(
            out_file,
            sep='\t',
            index=False,
            header=['wavelength','depth','depth_err'],
            float_format='%.6e',
            na_rep=''
        )
        print(f"✅ Combined spectrum written to {out_file}")

    def plot_overlay_with_errors(self):
        """
        Plot all spectra overlaid with 1-sigma error bars and save the figure.
        """
        plt.figure(figsize=(8,5))
        plt.errorbar(
            self.df_niriss['wavelength'], self.df_niriss['depth'],
            xerr=self.df_niriss['wavelength_err'], yerr=self.df_niriss['depth_err'],
            fmt='o-', label='NIRISS SOSS', color='C0', capsize=2
        )
        plt.errorbar(
            self.df_comb['wavelength'], self.df_comb['depth'],
            xerr=self.df_comb['wavelength_err'], yerr=self.df_comb['depth_err'],
            fmt='s--', label='Combined Model', color='C1', capsize=2
        )
        plt.errorbar(
            self.df_arch['wavelength'], self.df_arch['depth'],
            xerr=self.df_arch['wavelength_err'], yerr=self.df_arch['depth_err'],
            fmt='^:', label='Archival HST/Spitzer', color='C2', capsize=2
        )
        plt.errorbar(
            self.df_prism['wavelength'], self.df_prism['depth'],
            xerr=self.df_prism['wavelength_err'], yerr=self.df_prism['depth_err'],
            fmt='D-.', label='NIRSpec PRISM', color='C3', capsize=2
        )
        plt.xlabel('Wavelength (μm)')
        plt.ylabel('Transit Depth $(Rp/Rs)^2$')
        plt.title('Overlay of Spectra with $1\sigma$ Error Bars')
        plt.legend(loc='best', fontsize='small')
        plt.grid(alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(self.plots_dir, 'overlay_spectra_with_errors.png'), dpi=200)
        plt.close()
        print(f"✅ Saved plot with errors to {os.path.join(self.plots_dir, 'overlay_spectra_with_errors.png')}")

    def plot_overlay_no_errors(self):
        """
        Plot all spectra overlaid without error bars and save the figure.
        """
        plt.figure(figsize=(8,5))
        plt.plot(self.df_niriss['wavelength'], self.df_niriss['depth'], '-o', label='NIRISS SOSS', color='C0')
        plt.plot(self.df_comb['wavelength'], self.df_comb['depth'], '--s', label='Combined Model', color='C1')
        plt.plot(self.df_arch['wavelength'], self.df_arch['depth'], '-.^', label='Archival HST/Spitzer', color='C2')
        plt.plot(self.df_prism['wavelength'], self.df_prism['depth'], '-.D', label='NIRSpec PRISM', color='C3')
        plt.xlabel('Wavelength (μm)')
        plt.ylabel('Transit Depth $(Rp/Rs)^2$')
        plt.title('Overlay of Spectra (No Error Bars)')
        plt.legend(loc='best', fontsize='small')
        plt.grid(alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(self.plots_dir, 'overlay_spectra_no_errors.png'), dpi=200)
        plt.close()
        print(f"✅ Saved plot without errors to {os.path.join(self.plots_dir, 'overlay_spectra_no_errors.png')}")

    def plot_combined_spectrum(self):
        """
        Plot the final merged spectrum including error bars where available and save the figure.
        """
        plt.figure(figsize=(8,5))
        with_err = self.combined[self.combined['depth_err'].notna()]
        no_err = self.combined[self.combined['depth_err'].isna()]

        if not with_err.empty:
            plt.errorbar(
                with_err['wavelength'], with_err['depth'], yerr=with_err['depth_err'],
                fmt='o', label='With error bars', color='black', capsize=2
            )
        if not no_err.empty:
            plt.plot(
                no_err['wavelength'], no_err['depth'], 's',
                label='No error bar', color='blue'
            )
        plt.xlabel('Wavelength (μm)')
        plt.ylabel('Transit Depth $(Rp/Rs)^2$')
        plt.title('Combined Transmission Spectrum')
        plt.legend(loc='best', fontsize='small')
        plt.grid(alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(self.plots_dir, 'combined_spectrum.png'), dpi=200)
        plt.close()
        print(f"✅ Saved combined spectrum plot to {os.path.join(self.plots_dir, 'combined_spectrum.png')}")

    def run_all(self):
        """
        Execute all steps of the processing pipeline:
        load spectra, align baselines, merge spectra, and generate plots.
        """
        self.load_spectra()
        self.align_baselines()
        self.merge_spectra()
        self.plot_overlay_with_errors()
        self.plot_overlay_no_errors()
        self.plot_combined_spectrum()
