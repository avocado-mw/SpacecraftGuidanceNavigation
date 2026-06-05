# SpacecraftGuidanceNavigation
MAE182 Project

Included files
--------------
Files:
project_part_1.m     main script
project_dyn.m        dynamics + STM propagation
computed_A.m         9x9 linearized dynamics matrix
computed_Htilde.m    observation-state matrix + predicted observables
computed_obs.m       predicted range / range-rate only
project_obs_data.txt observation data file
project_setup.m
    Shared constants, initial conditions, covariance settings, and observation data.
project_batch.m
    Batch least-squares scaffold.
project_sequential.m
    Conventional sequential Kalman scaffold.
project_joseph.m
    Sequential Kalman scaffold with Joseph covariance update.
project_ekf.m
    Extended Kalman filter scaffold.
project_potter.m
    Potter square-root filter scaffold.
propagate_full_arc.m
propagate_step.m
plotEllipsoid.m
    Utility functions provided for organization and plotting.
project_obs_data.txt
    Observation file used throughout the project.
