%-------------------------------------------------------------------------%
%    project_ekf
%
%  DESCRIPTION:
%     Learner template for the extended Kalman filter (EKF).
%
%     In the EKF, the reference trajectory is updated after the measurement
%     update to reflect the current best estimate of the true trajectory:
%
%         X*_k <- X*_k + xhat_k
%
%     before propagating forward to the next observation time.
%
%  MAE 182 -- Spacecraft Guidance and Navigation
%-------------------------------------------------------------------------%
clearvars; close all; clc

%% 0) Shared setup
S = project_setup();

%% 1) TODO: implement the EKF
% Recommended strategy:
%   (a) begin from your sequential or Joseph script,
%   (b) after the measurement update, modify the reference state
%           Xref_k <- Xref_k + xhat_k
%   (c) reinitialize the next propagation with the updated reference state
%       and Phi(t_k,t_k) = I,
%   (d) decide whether you want to update immediately or after a warmup
%       interval once the filter has stabilized.
%
% Compare the EKF results against the batch and conventional sequential
% solutions.
