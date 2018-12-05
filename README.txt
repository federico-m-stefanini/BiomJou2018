README.txt

Source code for manuscript "A probabilistic network for the diagnosis of cardiopulmonary diseases"
by Alessandro Magrini, Davide Luciani and Federico M. Stefanini

For questions, comments or remarks about the code, please contact Alessandro Magrini (alessandro.magrini at outlook.com).

The codes have been written using the R software (platform: Windows 10.1, 64-bit).
R can be dowloaded from: https://cran.r-project.org/

To run MCMC simulation, analyze its output and implement the probabilistic network in GeNIe,
follow the instruction in the file 'CardiopulmNet_mcmc_START.r'.
The following R packages are required: 'coda', 'rjags', 'TeachingDemos'.
Software JAGS needs to be installed before loading the rjags package. It can be downloaded
from the official repository at: https://sourceforge.net/projects/mcmc-jags/

MCMC simulation output, figures and summaries will be saved in the folder 'save'.
Please make sure DO NOT change either its name, or the name of its subfolder.

The code to compute C-index values and the diagnosis of fictitious patient cases shown in
Tables 6 and 7 is not provided, because it was originally written for a previous open and
no more downloadable version of the software (2.0), which is not compatible with the current
one (2.1). At the moment, tables can be reproduced only using the graphical interface of
GeNIe 2.1, which can be downloaded from: https://www.bayesfusion.com/
