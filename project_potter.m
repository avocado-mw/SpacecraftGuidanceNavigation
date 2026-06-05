%-------------------------------------------------------------------------%
%    project_potter
%
%  DESCRIPTION:
%     Learner template for the Potter square-root update filter.
%
%     The Potter algorithm processes scalar observations one at a time and
%     updates the square-root factor W such that
%
%         P = W * W'
%
%     instead of updating P directly.
%
%  MAE 182 -- Spacecraft Guidance and Navigation
%-------------------------------------------------------------------------%
clearvars; close all; clc

%% 0) Shared setup
S = project_setup();

%% 1) TODO: initialize the square-root factor
% Wbar = chol( S.P0_bar );
%
% Then process the range and range-rate observations sequentially using the
% Potter scalar update equations from the notes.
%
% Suggested workflow at each observation epoch:
%   (a) time update the reference trajectory and Phi
%   (b) map the prior state deviation and square-root factor forward
%   (c) process the range observation
%   (d) process the range-rate observation
%   (e) recover P = W * W' if you wish to monitor the covariance trace
%
% Compare the Potter results against the conventional sequential filter and
% the Joseph formulation.
