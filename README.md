# acoustic-wave-propagation-model-

In this repository I share my code that solves the wave equation velocity-stress formulation in 2D acoustic media, using absorbing C-PML boundaries.

This serves to model the acoustic wave propagation phenomenon highly used in subsurface exploration processes.

Fourth order finite differences method were used in space and second order finite differences were used in time. 

This code also includes the implementation of pseudospectral method, which improves the results with finite differences.

The main file to run is "indata_Acustic_Marmousi.m"

In the file "a_forward_one_source.m" you can configure the geometry of sources and receivers, dx, dt, source wavelet, etc.


