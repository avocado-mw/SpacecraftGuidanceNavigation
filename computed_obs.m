%-------------------------------------------------------------------------%
%    computed_obs
%
%  DESCRIPTION:
%     Computes the predicted range and range-rate for one tracking station.
%
%     Inputs:
%       t   - time past epoch [s]
%       x   - full stored state [18x1]
%             [r; v; mu; J2; CD; station coordinates]
%       i   - station index in {1,2,3}
%
%     Outputs:
%       rho     - predicted range [m]
%       rhodot  - predicted range-rate [m/s]
%
%  NOTES:
%     - The satellite state is expressed in the inertial frame.
%     - Station coordinates are stored in ECEF and rotated into the inertial
%       frame through theta = thetadot * t.
%
%  MAE 182 -- Spacecraft Guidance and Navigation
%-------------------------------------------------------------------------%
function [ rho , rhodot ] = computed_obs( t , x , i )

%% 1) Earth rotation angle
% Relation between inertial and Earth-fixed frames.
thetadot = 7.2921158553e-5;    % [rad/s]
theta    = thetadot * t;
c        = cos( theta );
s        = sin( theta );

%% 2) Extract the selected station coordinates from the state vector
n  = 10 + 3 * ( i - 1 );
xi = x(n);
yi = x(n + 1);
zi = x(n + 2);

%% 3) Predicted range
rho = sqrt( ( x(1) - ( xi * c - yi * s ) )^2 + ...
            ( x(2) - ( xi * s + yi * c ) )^2 + ...
            ( x(3) - zi )^2 );

%% 4) Predicted range-rate
rhodot = 1 / rho * ( x(4) * ( x(1) - xi * c + yi * s ) + ...
                     x(5) * ( x(2) - yi * c - xi * s ) + ...
                     x(6) * ( x(3) - zi ) + ...
                     thetadot * ( ( x(1) * xi + x(2) * yi ) * s + ...
                                  ( x(1) * yi - x(2) * xi ) * c ) );

end
