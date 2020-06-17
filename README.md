# **RRL completeness utils**

This is an utility library to read, plot and recompute the RR Lyrae completeness maps presented in `Mateu et al. 2020 <>` (M20).

Gaia DR2 (VC+SOS), Pan-STARRS-1 and ASAS-SN-II completeness maps presented in M20 are available in csv format in the *maps* directory.

![see plot here](maps/VCSOS_final_completeness.png?raw=true "Gaia DR2 (VC+SOS) completeness map")


The Jupyter notebook *map_computing_examples.py* provides examples on how to compute 2D maps and completeness in the line of sight for a given catalogue and *map_plotting_examples.py* provide examples on how to plot the different precomputed maps published in `Mateu et al. 2020 <>`.

**REQUIREMENTS**

- Python modules required are ASTROPY, HEALPY, NUMPY, SCIPY, MATPLOTLIB is needed for plotting utilities.
