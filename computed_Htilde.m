%-------------------------------------------------------------------------%
%    computed_Htilde
%
%  DESCRIPTION:
%     Computes the observation-state matrix Htilde together with the
%     predicted range and range-rate for one tracking station.
%
%     Inputs:
%       t   - time past epoch [s]
%       x   - full stored state [18x1]
%             [r; v; mu; J2; CD; station coordinates]
%       i   - station index in {1,2,3}
%
%     Outputs:
%       Htilde  - 2x18 observation-state matrix
%       rho     - predicted range [m]
%       rhodot  - predicted range-rate [m/s]
%
%  NOTES:
%     Only the state components that appear explicitly in the observation
%     model have nonzero entries.  The partials with respect to mu, J2, and
%     CD are zero here because the observation equations depend directly only
%     on position, velocity, and station coordinates.
%
%  MAE 182 -- Spacecraft Guidance and Navigation
%-------------------------------------------------------------------------%
function [ Htilde , rho , rhodot ] = computed_Htilde( t , x , i )

%% 1) Earth rotation angle
thetadot = 7.2921158553e-5;    % [rad/s]
theta    = thetadot * t;
c        = cos( theta );
s        = sin( theta );

%% 2) Extract station coordinates from the full state
n  = 10 + 3 * ( i - 1 );
xi = x(n);
yi = x(n + 1);
zi = x(n + 2);

%% 3) Compute the predicted observables
rho = sqrt( ( x(1) - ( xi * c - yi * s ) )^2 + ...
            ( x(2) - ( xi * s + yi * c ) )^2 + ...
            ( x(3) - zi )^2 );

rhodot = 1 / rho * ( x(4) * ( x(1) - xi * c + yi * s ) + ...
                     x(5) * ( x(2) - yi * c - xi * s ) + ...
                     x(6) * ( x(3) - zi ) + ...
                     thetadot * ( ( x(1) * xi + x(2) * yi ) * s + ...
                                  ( x(1) * yi - x(2) * xi ) * c ) );

%% 4) Assemble Htilde = d[rho; rhodot]/dX
Htilde = zeros( 2 , 18 );
rinv   = 1 / rho;
d      = rhodot / rho;

% Partial derivatives of the range.
Htilde(1,1)   = rinv * ( x(1) - xi * c + yi * s );
Htilde(1,2)   = rinv * ( x(2) - yi * c - xi * s );
Htilde(1,3)   = rinv * ( x(3) - zi );
Htilde(1,n)   = rinv * ( xi - x(1) * c - x(2) * s );
Htilde(1,n+1) = rinv * ( yi - x(2) * c + x(1) * s );
Htilde(1,n+2) = rinv * ( zi - x(3) );

% Partial derivatives of the range-rate.
Htilde(2,1)   = rinv * ( x(4) + thetadot * ( xi * s + yi * c ) ...
                       - d * ( x(1) - xi * c + yi * s ) );
Htilde(2,2)   = rinv * ( x(5) + thetadot * ( yi * s - xi * c ) ...
                       - d * ( x(2) - yi * c - xi * s ) );
Htilde(2,3)   = rinv * ( x(6) - d * ( x(3) - zi ) );
Htilde(2,4)   = rinv * ( x(1) - xi * c + yi * s );
Htilde(2,5)   = rinv * ( x(2) - yi * c - xi * s );
Htilde(2,6)   = rinv * ( x(3) - zi );
Htilde(2,n)   = rinv * ( -x(4) * c - x(5) * s ...
                       + thetadot * ( x(1) * s - x(2) * c ) ...
                       - d * ( xi - x(1) * c - x(2) * s ) );
Htilde(2,n+1) = rinv * ( -x(5) * c + x(4) * s ...
                       + thetadot * ( x(2) * s + x(1) * c ) ...
                       - d * ( yi - x(2) * c + x(1) * s ) );
Htilde(2,n+2) = rinv * ( -x(6) - d * ( zi - x(3) ) );

end
