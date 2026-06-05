%-------------------------------------------------------------------------%
%    sequential_gap_reset
%
%  DESCRIPTION:
%     Reset the sequential deviation state and covariance at the start of a
%     new visibility pass when the elapsed time since the last observation
%     exceeds the project gap threshold.
%
%  MAE 182 -- Spacecraft Guidance and Navigation
%-------------------------------------------------------------------------%
function [ x_hat_prev , P_prev ] = sequential_gap_reset( S , dt , x_hat_prev , P_prev )

if dt > S.gapResetTime
    x_hat_prev = S.x0_bar;
    P_prev     = S.P0_bar;
end

end
