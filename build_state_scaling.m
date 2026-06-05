%-------------------------------------------------------------------------%
%    build_state_scaling
%
%  DESCRIPTION:
%     Diagonal scaling matrix S for state normalization.  Physical states
%     relate to normalized coordinates through  x = S * z .
%
%     Scales are chosen so all components of H*z are O(1):
%       r [m]           -> 1e6
%       v [m/s]         -> 1e3
%       mu [m^3/s^2]    -> 1e14
%       J2 [-]          -> 1e-3
%       CD [-]          -> 1
%       station coords  -> 1e6
%
%  MAE 182 -- Spacecraft Guidance and Navigation
%-------------------------------------------------------------------------%
function S = build_state_scaling()

scales = [ 1e6 , 1e6 , 1e6 , ...   % r
           1e3 , 1e3 , 1e3 , ...   % v
           1e14 , ...              % mu
           1e-3 , ...              % J2
           1 , ...                 % CD
           1e6 , 1e6 , 1e6 , ...   % station 101
           1e6 , 1e6 , 1e6 , ...   % station 337
           1e6 , 1e6 , 1e6 ];      % station 394

S = diag( scales );

end
