# Subevents2019


Scripts that were used in the following publication:
Danre, Yin, Lipovsky, Denolle "Earthquakes within earthquakes: Patterns in rupture complexity" (in revision in GRL 2019).

This directory contains python scripts to perform a subevent decomposition of source time functions using Gaussian kernels. It contains resampled USGS data set of 180 STFs (as published in Hayes 2017), 500 STFs simulated using SBIEM from Pablo Ampuero.
The directory also contains matlab scripts used to run SBIEM in the context of the publication.


# requirements.
numpy,matplotlib,sys,os,glob,scipy

# List of scripts

subevent_scardec.py: main script for Danre et al, 2019. Includes plotting for the figures in the paper.
test_triangle.py:       copy of main script but with triangle kernels
test_usgs.py:           same as in main but for USGS database.
test_simulations.py     same as in main for for simulated STFs.

# list of STF files:
USGS_STF/ * . newstf:  all USGS STFs sampled at dt=0.0730s, 180 of them (1990-2017)
SIM_STF/ * .dat     :  all simulated STFs sampled at dt=0.0730s


# list of variables
allvar.npz:   results out of main
allvar_triangle.npz:  results out of test_triangle.py
allvar_usgs.npz:    results out of test_usgs.py
allvar_sim.npz:     results out of test_simulations.py



Contact marine for any question: mdenolle(AT)fas(DOT)harvard(DOT)edu.
