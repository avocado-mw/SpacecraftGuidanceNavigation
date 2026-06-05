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
project_obs_data.txt observation data file (all three stations)
project_obs_data_337.txt
    Station 337 only; used conceptually for sequential filters.
project_setup.m
    Shared constants, initial conditions, covariance settings, and observation data.
project_batch.m
    Batch least-squares processor (all stations).
project_sequential.m
    Conventional sequential Kalman processor (station 337 only).
project_joseph.m
    Sequential Kalman processor with Joseph covariance update.
project_ekf.m
    Extended Kalman filter with reference-trajectory updates.
project_potter.m
    Potter square-root sequential filter.
potter_scalar_update.m
    Scalar Potter measurement update helper.
plot_filter_residuals.m
    Residual plotting versus epoch time [s].
sequential_gap_reset.m
    Resets covariance at visibility-pass boundaries.
propagate_full_arc.m
propagate_step.m
plotEllipsoid.m
    Utility functions provided for organization and plotting.

Observation data notes
----------------------
The full file interleaves three ground stations (101, 337, 394) with multi-hour
visibility gaps. Sequential filters use station 337 only via
`project_setup('sequential')` because mixing stations in one sequential pass
causes unstable covariance growth and obscures residual trends. Batch least
squares continues to use all stations.

Running the project
-------------------
In MATLAB or Octave, run the scripts from this directory:

    project_part_1      % Part I: propagation, residuals, RMS
    project_batch       % batch least squares (3 iterations, all stations)
    project_sequential  % conventional sequential filter (station 337)
    project_joseph      % Joseph-form sequential filter
    project_ekf         % extended Kalman filter
    project_potter      % Potter square-root filter

Each estimator saves a `.mat` results file and plots residuals versus epoch
time when `S.makePlots` is true in `project_setup.m`.
