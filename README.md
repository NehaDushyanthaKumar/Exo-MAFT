# Exo-MAFT

[![License: MIT](https://cdn.prod.website-files.com/5e0f1144930a8bc8aace526c/65dd9eb5aaca434fac4f1c34_License-MIT-blue.svg)](/LICENSE) [![Python 3.10](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/downloads/release/python-360/) [![A rectangular badge, half black half purple containing the text made at Code Astro](https://img.shields.io/badge/Made%20at-Code/Astro-blueviolet.svg)](https://semaphorep.github.io/codeastro/)

The Multi-Wavelength Atmospheric Feature Tracker for Exoplanets (Exo-MAFT) addresses a critical limitation in exoplanet atmospheric studies: existing tools predominantly focus on single-wavelength analysis, whereas atmospheric features vary significantly across wavelengths and require coordinated multi-spectral analysis. This package will help researchers to track atmospheric features over multiple wavelength ranges, identify correlations across spectra, and optimize observation strategies.

![IMG_0645](https://github.com/user-attachments/assets/7ed7116e-55a1-4f41-969f-0d8d5ed8e7a4)

## Installation

You can install **Exo-MAFT** using pip:

```bash
pip install exomaft

```
We prefer for the users to git clone instead (to download the large linelist files):

```bash
git clone https://github.com/NehaDushyanthaKumar/Exo-MAFT.git
```

## Quick Start – Example: WASP-39b

You're now ready to use the **Exo-MAFT** package! Below is an example workflow for **WASP-39b** using publicly available JWST data.  
You can adapt the file paths and inputs for **other exoplanets** that have been observed through **transmission spectroscopy**.

### In Python

```python
import exomaft as maft

```
### Change to your file locations
```
file_paths = {
    'niriss': 'exomaft/data/NIRISS SOSS/FIREFly_transmission_spectrum_R_300_order_1.txt',
    'combined': 'exomaft/data/Combined/BestFitModel.txt',
    'archival': 'exomaft/data/Combined/Archival_Data.txt',
    'prism': 'exomaft/data/NIRSpec_PRISM/FIREFLy_transit_spectrum.txt'
}
```


### Change to your file locations
```
file_paths = {
    'niriss': 'exomaft/data/NIRISS SOSS/FIREFly_transmission_spectrum_R_300_order_1.txt',
    'combined': 'exomaft/data/Combined/BestFitModel.txt',
    'archival': 'exomaft/data/Combined/Archival_Data.txt',
    'prism': 'exomaft/data/NIRSpec_PRISM/FIREFLy_transit_spectrum.txt'
}
```

### Process and combine transmission spectra from multiple instruments
```
processor = maft.TransmissionSpectrumProcessor(file_paths=file_paths)
processor.run_all()
```

### Identify molecular feature lines
```
tracker = maft.FeatureTracker()
tracker.plot_and_save()
```

### Run cross-correlation analysis
```
correlator = maft.CrossCorrelator()
correlator.run_cross_correlation()
```







