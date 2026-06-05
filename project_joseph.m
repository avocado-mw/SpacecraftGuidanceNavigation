%-------------------------------------------------------------------------%
%    project_joseph
%
%  DESCRIPTION:
%     Learner template for the sequential processor using the Joseph-form
%     covariance measurement update.
%
%     This is almost identical to project_sequential.m, except that the
%     covariance update is
%
%         P_k = (I-KH) Pbar_k (I-KH)' + K R K'
%
%     which is numerically more robust.
%
%  MAE 182 -- Spacecraft Guidance and Navigation
%-------------------------------------------------------------------------%
clearvars; close all; clc

%% 0) Shared setup
S = project_setup();

%% 1) TODO: adapt your conventional sequential script
% Start from project_sequential.m and replace only the covariance update
% with the Joseph form.  Compare the residuals and covariance traces.
