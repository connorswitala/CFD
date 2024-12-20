# CFD
Computational Fluid Dynamics code for both inviscid and viscous.

% 'Flat_Plate_Generator.m', 'Cylinder_Grid_Generator.m', 'Ramp_Grid_Generator.m' creates the computational domain and outputs geometrical quantities of face normals, side areas, volume, number of cells in x direction, number of cells in y direction, cell vertices, cell centers. 

% 'Inviscid_Ramp_DPLR.m' and 'Inviscid_Cylinder_DPLR.m' are the main files that run the grid generators, and iterate using DPLR until convergence for inviscid flow. 

% 'Viscous_DPLR.m' includs viscous effects.

% The code is not completely modularized yet so each case needs to be run separately. It is also extremely inefficent at the moment and totally unreadable. 


