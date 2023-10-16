# Code for "Energy coupling and stoichiometry of Zn<sup>2+</sup>/H<sup>+</sup> antiport by the cation diffusion facilitator YiiP"

[![DOI](https://zenodo.org/badge/383897104.svg)](https://doi.org/10.5281/zenodo.8357619)



This repository contains code used for some of the computational aspects of the paper

* Adel Hussein, Shujie Fan, Maria Lopez-Redondo, Ian Kenney, Xihui Zhang, Oliver Beckstein, and David L. Stokes. Energy coupling and stoichiometry of Zn<sup>2+</sup>/H<sup>+</sup> antiport by the cation diffusion facilitator YiiP. eLife, Apr 2023. doi: [10.7554/elife.87167](https://doi.org/10.7554/elife.87167)

Code and data are archived: 
* Raw data files are archived in OSF project https://osf.io/r95qu (DOI [10.17605/OSF.IO/Y8BA2](https://doi.org/10.17605/OSF.IO/Y8BA2))
* This GitHub code repository is archived in Zenodo under DOI [10.5281/zenodo.8357619](https://doi.org/10.5281/zenodo.8357619)

## The following packages are required to use these scripts.
* [NumPy](https://numpy.org/)
* [SciPy](https://scipy.org/)
* [Matplotlib](https://matplotlib.org/)
* [MDAnalysis](https://www.mdanalysis.org/)
* [tqdm](https://github.com/tqdm/tqdm)
* [seaborn](https://seaborn.pydata.org/)
* [pandas](https://pandas.pydata.org/)
* [multibind](https://github.com/Becksteinlab/multibind)
* [CpHMD-Analysis](https://gitlab.com/shenlab-amber-cphmd/cphmd-analysis) (see https://github.com/Hendejac/CpHMD-Analysis for the original `cphmdanalysis.py` script/module)

## analysis
`angle.py`, `rms.py`, and `R210_D72.py` are python scripts used to analyze equilibrium simulations, taking `md.gro` and `md[].xtc` as inputs, while the reference for this analysis is derived from the Protein Data Bank, accessible at https://www.rcsb.org/structure/5VRF.

`cphmd.py` is the python script for analyzing CpHMD simulations, with `.log` and `.lamb` files as inputs.

## MST inference
The folders `siteA` and `siteB` each contain scripts and corresponding example input data for performing MST inference on site A and site B, respectively. While the scripts themselves are identical for both site A and site B, the distinction lies in their input and output files.

`monte_carlo.py` serves as the script responsible for executing a Monte Carlo process to refine the pKa values and binding free energies. It requires a `.csv` file that defines the states within the thermodynamic model, exemplified by `state.csv`, and an initial file, following the same format as `result.csv`. The `.csv` files located in the `target` folder as input data serves as the reference values for the Monte Carlo scheme. These files are reformatted experimental MST data, representing the data in Figure 4 of the paper.

`analysis.ipynb` is a jupyter-notebook showing how to use the output of the monte_carlo process, `result.csv`, to generate the figures in `Figure 4 - figure supplement 5`, and `Figure 4 - figure supplement 6`.The `result.csv` files in `siteA` and `siteB` folders are reformatted versions of the supplementary data `figure4_figure-supplement5_source-data3.csv` and `figure4_figure-supplement6_source-data3.csv`, respectively.
