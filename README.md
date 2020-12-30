
# Chance-Constrained Programming for Autonomous Vehicles in Uncertain Environments
This repository contains the source code used for the Semester Thesis titled "Chance-Constrained Programming for Autonomous Vehicles in Uncertain Environments" submitted at ETH Zurich under the supervision of prof. Maryam Kamgarpour. DOI: [10.3929/ethz-b-000272614](https://doi.org/10.3929/ethz-b-000272614)

## Dependencies
The following dependencies need to be installed/configured and must be on the MATLAB path:
- [YALMIP](https://yalmip.github.io/) (configured with an appropriate solver like [CPLEX]( https://www.ibm.com/analytics/cplex-optimizer))
- [matlab2tikz](https://github.com/matlab2tikz/matlab2tikz)
- `utils` folder

## Generating plots
In order to reproduce the simulations and plots of a specific chapter, navigate inside the corresponding chapter folder and run the `generatePlotsAll` MATLAB script. The script will run all the necessary simulations and produce all the plots of that chapter. The plots will also be saved under a `plots` folder in TikZ format, which can then be readily included in a LaTeX document.

## Citing this work
Please cite the original thesis when using any part of this code. BibTeX citation data:
```
@misc{Lefkopoulos2018,
	author	  	 = "V. Lefkopoulos and M. Kamgarpour",
	title	  	 = "Chance-Constrained Programming for Autonomous Vehicles in Uncertain Environments",
	howpublished     = "ETH Zurich", 
	year	  	 = "2018",
}
```
