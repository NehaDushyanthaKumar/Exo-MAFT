README for feature_tracking.py

Overview
feature_tracking.py is a standalone script that:
1. Reads a combined transmission spectrum ( data/combined_spectrum.txt ).
2. Reads a HITRAN‐style line list ( .par format) via fixed-width parsing.
3. Selects user-specified molecules (e.g. H₂O, CO₂, NH₃, CO, CH₄, HCN).
4. Clusters each molecule’s top‑N strongest transitions into wavelength bins.
5. Plots the spectrum with semi‑transparent vertical lines marking molecular features.
6. Saves the marker wavelengths into an output text file ( output/feature_lines.txt ).
7. Exports the figure to plots/feature_tracking.png .

Prerequisites
• Python 3.7+ with the following packages:
• numpy
• pandas
• matplotlib

Getting the Linelist
The script expects a HITRAN-style .par file containing molecular line data. You can generate or filter
your own line list using the getline.py utility in this repository:

python getline.py --molecule H2O --output data/H2O.par
Adjust --molecule to any species supported by HITRAN (e.g. CO2, NH3, CO, CH4, HCN).

Usage
1. Place your combined spectrum in data/combined_spectrum.txt .
2. Place or generate your line list in data/6894c8ca.par (or any filename you prefer).
3. Adjust the COMBINED_SPECTRUM and PAR_FILE constants at the top of
feature_tracking.py if needed.
4. Run the script:
python feature_tracking.py
This will create:
• plots/feature_tracking.png : a figure showing the spectrum and molecular feature lines.
• output/feature_lines.txt : a tab‑delimited table listing each molecule and its marker
wavelengths.

1

Script Breakdown
# CONFIGURATION
#
Defines input/output paths and which molecules (by HITRAN IDs) to track.
COMBINED_SPECTRUM = 'data/combined_spectrum.txt'
PAR_FILE
= 'data/6894c8ca.par'
MOLECULES = {'H2O':1, 'CO2':2, ...}
TOP_N
= 30
# how many strongest lines to cluster per molecule
BIN_WIDTH
= 0.02
# μm, cluster width
ALPHA_LINES = 0.3
# transparency for feature lines
# LOAD COMBINED SPECTRUM
#
Reads wavelength, depth, and depth_err into numpy arrays.
spec = pd.read_csv(...)
wave = spec['wavelength_um'].values
depth = spec['depth'].values
# LOAD LINELIST (.par)
#
Uses pandas.read_fwf with fixed-width columns defined in README
#
Casts to numeric and computes 'wavelength_um = 1e4/nu'.
def load_linelist_par(path):
colspecs = [...]
names = [...]
df = pd.read_fwf(path, colspecs=colspecs, names=names)
df['wavelength_um'] = 1e4 / df['nu']
return df
# PLOT SETUP
#
Opens a matplotlib figure, plots the spectrum and shaded error band.
fig, ax = plt.subplots()
ax.plot(wave, depth)
ax.fill_between(...)
# OVERLAY FEATURE LINES
#
For each molecule:
#
• select top TOP_N lines by strength 'sw'
#
• cluster them into BIN_WIDTH bins
#
• compute cluster centers
#
• draw semi-transparent vertical lines at those centers
#
• record (molecule, center) in an output list
# WRITE OUTPUT TABLE
#
Saves a tab-delimited file 'output/feature_lines.txt'
df_out = pd.DataFrame(out_rows, columns=['molecule','wavelength_um'])
df_out.to_csv(OUT_LINES, sep='\t', index=False)
# FINALIZE PLOT

2

#
Add labels, legend, grid, then save to 'plots/feature_tracking.png'.
plt.savefig(PLOT_FILE)
For detailed cross-correlation analysis, see cross_correlation.py and its accompanying README.
For any questions or customization, refer to the in-line code comments or open an issue on GitHub.

3

