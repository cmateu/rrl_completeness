#  RRL completeness utils

This is an utility library to read, plot and recompute the RR Lyrae completeness maps presented in [Mateu et al. 2020](https://arxiv.org/abs/2006.09416).

Gaia DR2 (VC+SOS), Pan-STARRS-1 and ASAS-SN-II completeness maps presented in M20 are available in csv format in the *maps* directory.

![see plot here](maps/VCSOS_final_completeness.png?raw=true "Gaia DR2 (VC+SOS) completeness map")

----------
### EXAMPLES

- Utility functions are available in the completeness_util module

- The Jupyter notebook  [map_plotting_examples.py](https://github.com/cmateu/rrl_completeness/blob/master/examples/map_plotting_examples.ipynb) shows how to plot the different precomputed maps published in [Mateu et al. 2020](https://arxiv.org/abs/2006.09416).
- The Jupyter notebook [map_computing_examples.py](https://github.com/cmateu/rrl_completeness/blob/master/examples/map_computing_examples.ipynb) shows how to compute 2D maps and completeness in the line of sight for a given catalogue.

----------

### REQUIREMENTS

- Python modules required are NUMPY, SCIPY, ASTROPY AND HEALPY. PANDAS, MATPLOTLIB and [optionally SEABORN] are needed for plotting utilities.

----------

### INSTALLATION

In a terminal, run the following command:

    sudo python setup.py install

and source your .cshrc

If you do not have root access, you can install in the custom directory path_to_dir.
First, add the directory's path path_to_dir and path_to_dir/lib/python2.7/site-packages/
to the PYTHONPATH variable in your .cshrc (or .bashrc) file and source it. Then install using the --prefix option::

    python setup.py install --prefix=path_to_dir

Add path_to_dir/bin to your PATH in your .csrhc or .bashrc file.
