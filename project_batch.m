%-------------------------------------------------------------------------%
%    project_batch
%
%  DESCRIPTION:
%     Learner template for the batch least-squares processor used in the
%     MAE 182 term project.
%
%     This script assumes that you have already implemented:
%         project_dyn.m
%         computed_A.m
%         computed_Htilde.m
%         computed_obs.m
%
%     The batch processor should:
%       (1) propagate the reference trajectory and STM over the full data arc
%       (2) build H_i = Htilde_i * Phi(t_i,t_0)
%       (3) accumulate Lambda and N
%       (4) solve the normal equations for x0_hat
%       (5) iterate the nominal state update several times
%       (6) report pre-fit / post-fit residuals and the position covariance
%
%  MAE 182 -- Spacecraft Guidance and Navigation
%-------------------------------------------------------------------------%
clearvars; close all; clc

%% 0) Shared setup
S = project_setup();

%% 1) Initialize the batch iteration
X0_ode = S.X0_ode;
x0_bar = S.x0_bar;
P0_bar = S.P0_bar;

nIter        = S.batch_nIter;
x0_hat_hist  = zeros( 18 , nIter );
rho_rms_hist = zeros( nIter , 1 );
rhodot_rms_hist = zeros( nIter , 1 );

for pass = 1:nIter

    %% 2) Propagate the reference trajectory and resampled STM history
    [ t20 , x20 , Phi20 ] = propagate_full_arc( S , X0_ode );

    %% 3) Initialize the normal equations
    Lambda = P0_bar \ eye( 18 );
    N      = Lambda * x0_bar;

    % Preallocate storage for this pass
    nObs   = size( S.Y , 1 );
    rho    = zeros( nObs , 1 );
    rhodot = zeros( nObs , 1 );
    y      = zeros( 2 , nObs );
    H      = zeros( 2 , 18 , nObs );

    %% 4) TODO: loop over the observations
    % For each observation:
    %   (a) recover the resampled index j = S.Yidx(k,1)
    %   (b) compute Htilde, rho, rhodot using computed_Htilde
    %   (c) form H(:,:,k) = Htilde * Phi20(:,:,j)
    %   (d) build the observation residual y(:,k) = Yobs - Ycomp
    %   (e) accumulate
    %
    %           Lambda = Lambda + H' * W * H
    %           N      = N      + H' * W * y
    %
    % Use the project handout and lecture notes as your guide.

    %% 5) TODO: compute RMS values for this pass
    % rho_rms_hist(pass)    = ...
    % rhodot_rms_hist(pass) = ...

    %% 6) TODO: solve the normal equations and update the nominal state
    % x0_hat = ...
    % P0     = ...
    %
    % Then update:
    %   X0_ode(1:18) = X0_ode(1:18) + x0_hat
    %   x0_bar       = x0_bar - x0_hat
    %
    % Store x0_hat in x0_hat_hist(:,pass).

end

%% 7) TODO: add the plots and summary outputs required by the project
% At a minimum, the final project should show:
%   - pre-fit and post-fit residuals
%   - RMS values
%   - the position-error ellipsoid
%   - discussion of the fixed station and covariance choices
