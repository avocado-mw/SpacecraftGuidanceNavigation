%-------------------------------------------------------------------------%
%    project_sequential
%
%  DESCRIPTION:
%     Learner template for the conventional sequential (Kalman) processor.
%
%     The key sequential steps are:
%       (1) propagate the reference trajectory and STM from t_{k-1} to t_k
%       (2) compute the time update for xbar_k and Pbar_k
%       (3) form the observation deviation y_k and Htilde_k
%       (4) compute the Kalman gain
%       (5) perform the measurement update for xhat_k and P_k
%
%  MAE 182 -- Spacecraft Guidance and Navigation
%-------------------------------------------------------------------------%
clearvars; close all; clc

%% 0) Shared setup
S = project_setup();

%% 1) Initial conditions for the sequential processor
x_hat_prev = S.x0_bar;      % state deviation estimate
P_prev     = S.P0_bar;
Xref0_ode  = S.X0_ode;
Phi0       = S.Phi0;

nObs = size( S.Y , 1 );

% Preallocate storage
x_hat   = zeros( 18 , nObs );
P_store = zeros( 18 , 18 , nObs );
pre_fit = zeros( 2 , nObs );
post_fit = zeros( 2 , nObs );

%% 2) TODO: process each observation in order
% For k = 1:nObs:
%   (a) propagate the reference trajectory and Phi from t_{k-1} to t_k
%   (b) time update:
%           x_bar = Phi_k * x_hat_prev
%           P_bar = Phi_k * P_prev * Phi_k'
%   (c) compute y_k = Y_k - G(X*_k, t_k) and Htilde_k
%   (d) Kalman gain:
%           K = P_bar * Htilde_k' / ( Htilde_k * P_bar * Htilde_k' + R )
%   (e) pre-fit residual:
%           pre_fit = y_k - Htilde_k * x_bar
%   (f) state update:
%           x_hat_k = x_bar + K * pre_fit
%   (g) covariance update (conventional form):
%           P_k = ( I - K * Htilde_k ) * P_bar
%   (h) post-fit residual:
%           post_fit = y_k - Htilde_k * x_hat_k
%   (i) set x_hat_prev = x_hat_k, P_prev = P_k, and
%       reinitialize Xref0_ode = [ X*_k ; Phi0 ] for the next step.
%
% IMPORTANT:
%   This script is the conventional sequential estimator.
%   The reference trajectory is NOT updated with x_hat_k here.

%% 3) TODO: compute RMS values and produce the requested plots
