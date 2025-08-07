README for cross_correlation.py

Overview
cross_correlation.py is a standalone script that:
1. Reads a combined transmission spectrum ( data/combined_spectrum.txt ).
2. Reads a HITRAN‐style line list ( .par format) via fixed-width parsing.
3. Selects user-specified molecules (e.g. H₂O, CO₂, NH₃, CO, CH₄, HCN).
4. Builds binary templates of each molecule’s strongest line centers.
5. Computes the cross-correlation function (CCF) between the observed spectrum and each
template.
6. Plots all CCF traces over wavelength to show where each molecule’s pattern best matches data.
7. Saves the peak CCF values into an output file ( output/ccf_peaks.txt ).
8. Exports the figure to plots/ccf_all.png .

Prerequisites
• Python 3.7+ with these packages:
• numpy
• pandas
• matplotlib
• scipy

Getting the Linelist
You can generate or filter your own HITRAN‐style .par line lists using the getline.py utility in this
repository:

python getline.py --molecule H2O --output data/H2O.par
Replace --molecule with any supported species (CO2, NH3, CO, CH4, HCN).

Usage
1. Place your combined spectrum in data/combined_spectrum.txt .
2. Generate or place your line list in data/6894c8ca.par (or update the PAR_FILE path).
3. Adjust the COMBINED_SPECTRUM and PAR_FILE constants at the top of
cross_correlation.py if needed.
4. Run:
python cross_correlation.py
This will create:
• plots/ccf_all.png : overlaid CCF traces (one per molecule) vs. wavelength.
• output/ccf_peaks.txt : tab-delimited table of each molecule and its peak CCF value.

1

Script Breakdown
# CONFIGURATION
#
Defines input/output paths and which molecules to track by HITRAN ID.
COMBINED_SPECTRUM = 'data/combined_spectrum.txt'
PAR_FILE
= 'data/6894c8ca.par'
MOLECULES = {'H2O':1, 'CO2':2, 'NH3':3, 'CO':5, 'CH4':6, 'HCN':8}
TOP_N
= 30
# number of top lines per molecule
BIN_WIDTH = 0.02

# μm clustering width for line centers

# MAKE OUTPUT DIRECTORIES
os.makedirs('plots', exist_ok=True)
os.makedirs('output', exist_ok=True)
# LOAD SPECTRUM
#
Reads wavelength, depth, depth_err into numpy arrays for CCF.
spec = pd.read_csv(...)
wave = spec['wavelength_um'].values
depth = spec['depth'].values
# LOAD LINELIST (.par)
#
Parses fixed-width columns, casts to numeric, computes wavelength.
def load_linelist_par(path):
colspecs = [...] # fixed-width ranges from readme-1.txt
names
= [...] # column names mapping
df = pd.read_fwf(path, colspecs=colspecs, names=names)
df['wavelength_um'] = 1e4 / df['nu']
return df
lines = load_linelist_par(PAR_FILE)
# UTILITY: z-score arrays for correlation
def zscore(arr):
return (arr - arr.mean()) / arr.std()
# CROSS-CORRELATE & RECORD PEAKS
results = []
for name, mol_id in MOLECULES.items():
dfm = lines[lines['molec_id']==mol_id]
top_waves = dfm.nlargest(TOP_N, 'sw')['wavelength_um'].values
# cluster into bins and compute centers
centers = cluster_centers(top_waves, BIN_WIDTH)
# build template spikes, zscore, and compute ccf
template = build_binary_template(wave, centers)
ccf = correlate(zscore(depth), zscore(template), mode='same')
peak = ccf.max()
results.append((name, peak))
plt.plot(wave, ccf, label=f"{name} (peak {peak:.1f})")
# SAVE PEAKS TO TEXT

2

pd.DataFrame(results, columns=['molecule','ccf_peak']).to_csv(
'output/ccf_peaks.txt', sep='\t', index=False, float_format='%.3f'
)
# FINALIZE PLOT
plt.xlabel('Wavelength (μm)')
plt.ylabel('CCF')
plt.title('Cross-Correlation Functions')
plt.legend()
plt.savefig('plots/ccf_all.png')
plt.show()

Interpreting Outputs
Numeric Peaks ( ccf_peaks.txt )
• Each ccf_peak is the maximum correlation between your spectrum and that molecule’s
template.
• Higher peaks indicate stronger matches (i.e., that molecule’s spectral lines are more clearly
present).

CCF Plot ( ccf_all.png )
• X-axis: Wavelength (µm).
• Y-axis: Correlation strength (z-scored).
• Traces: One colored line per molecule shows how well its template aligns at each wavelength.
• Legend: Molecule names with their peak CCF values.
Use these outputs to rank molecular detections, identify specific band locations, and guide follow-up
analysis.

3

