# SpacecraftGuidanceNavigation
MAE182 Project

Included files
--------------
Files:
project_part_1.m     main script (Part I reference solution)
project_dyn.m        dynamics + STM propagation
computed_A.m         9x9 linearized dynamics matrix
computed_Htilde.m    observation-state matrix + predicted observables
computed_obs.m       predicted range / range-rate only
project_obs_data.txt observation data file
project_setup.m
    Shared constants, initial conditions, covariance settings, and observation data.
project_batch.m
    Batch least-squares processor.
project_sequential.m
    Conventional sequential Kalman processor.
project_joseph.m
    Sequential Kalman processor with Joseph covariance update.
project_ekf.m
    Extended Kalman filter with reference-trajectory updates.
project_potter.m
    Potter square-root sequential filter.
potter_scalar_update.m
    Scalar Potter measurement update helper.
propagate_full_arc.m
propagate_step.m
plotEllipsoid.m
    Utility functions provided for organization and plotting.
project_obs_data.txt
    Observation file used throughout the project.

Running the project
-------------------
In MATLAB or Octave, run the scripts from this directory:

    project_part_1      % Part I: propagation, residuals, RMS
    project_batch       % batch least squares (3 iterations)
    project_sequential  % conventional sequential filter
    project_joseph      % Joseph-form sequential filter
    project_ekf         % extended Kalman filter
    project_potter      % Potter square-root filter

Each estimator saves a `.mat` results file and produces residual / ellipsoid plots when `S.makePlots` is true in `project_setup.m`.
