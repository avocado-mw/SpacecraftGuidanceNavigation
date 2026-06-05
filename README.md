# SpacecraftGuidanceNavigation
MAE182 Project

## Running the project

In MATLAB or Octave, run from this directory:

    project_part_1
    project_batch
    project_sequential
    project_joseph
    project_ekf
    project_potter
    project_sensitivity_p0   % P0 station-variance sensitivity study

## State and covariance setup

The 18-state estimation vector is

    [ r(3); v(3); mu; J2; CD; x_s1(3); x_s2(3); x_s3(3) ]

Station 101 (`x_s1`) is treated as known (`P0 = 1e-10`). Stations 337 and 394
(`x_s2`, `x_s3`) are unknown with `P0 = 1e6` as specified in the handout.

Measurement noise is

    R = diag( (0.01 m)^2, (0.001 m/s)^2 )

Batch and sensitivity processors normalize states as `x = S*z` and solve the
stacked least-squares system with QR decomposition instead of `inv(H'WH)`.
Sequential filters in the sensitivity study use the same normalization.
Joseph-form covariance updates are used by default.
