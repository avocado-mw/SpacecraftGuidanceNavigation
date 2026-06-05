%-------------------------------------------------------------------------%
%    computed_A
%
%  DESCRIPTION:
%     Computes the linearized dynamics matrix A for the 9-state propagated
%     subset used in MAE 182 Project -- Part I:
%
%         x_dyn = [ r ; v ; mu ; J2 ; CD ]
%
%     The model contains two-body gravity, J2 oblateness, and atmospheric
%     drag.  Rows/columns 7--9 correspond to the estimated parameters
%     [mu, J2, CD].  Since these parameters are modeled as constants, the
%     last three rows are zero.
%
%  INPUT:
%     x = full stored state [18x1] or longer, whose first 9 entries are
%         [r; v; mu; J2; CD]
%
%  OUTPUT:
%     A = 9x9 Jacobian matrix dF/dx_dyn
%
%  MAE 182 -- Spacecraft Guidance and Navigation
%-------------------------------------------------------------------------%
function A = computed_A( x )

%% 1) Constants and estimated parameters
Re = 6378136.3;        % Earth radius [m]
mu = x(7);             % Earth gravitational parameter [m^3/s^2]
J2 = x(8);             % Earth oblateness coefficient [-]
Cd = x(9);             % Drag coefficient [-]

S        = 3.0;        % spacecraft projected area [m^2]
m        = 970;        % spacecraft mass [kg]
rho0     = 3.614e-13;  % reference atmospheric density [kg/m^3]
r0       = 7e5 + Re;   % reference radius for density model [m]
H        = 88667.0;    % scale height [m]
thetadot = 7.2921158553e-5;  % Earth rotation rate [rad/s]

%% 2) State aliases for readability
x1 = x(1);  x2 = x(2);  x3 = x(3);
x4 = x(4);  x5 = x(5);  x6 = x(6);

r    = norm( x(1:3) );
rhoA = rho0 * exp( -( r - r0 ) / H );
Va   = [ x4 + thetadot * x2 ; ...
         x5 - thetadot * x1 ; ...
         x6 ];
VA   = norm( Va );

%% 3) Allocate Jacobian
A      = zeros( 9 , 9 );
A(1,4) = 1;
A(2,5) = 1;
A(3,6) = 1;

%% 4) Position-acceleration and parameter-acceleration partials
A(4,1) = -( mu / r^3 ) * ( 1 - (3/2) * J2 * ( Re / r )^2 * ( 5 * ( x3 / r )^2 - 1 ) ) + ...
          3 * ( mu * x1^2 / r^5 ) * ( 1 - (5/2) * J2 * ( Re / r )^2 * ( 7 * ( x3 / r )^2 - 1 ) ) + ...
          (1/2) * Cd * ( S / m ) * rhoA * VA * ( x1 / ( r * H ) ) * ( x4 + thetadot * x2 ) - ...
          (1/2) * Cd * ( S / m ) * ( rhoA / VA ) * ( -thetadot * x5 + thetadot^2 * x1 ) * ( x4 + thetadot * x2 );

A(4,2) = 3 * ( mu * x1 * x2 / r^5 ) * ( 1 - (5/2) * J2 * ( Re / r )^2 * ( 7 * ( x3 / r )^2 - 1 ) ) + ...
          (1/2) * Cd * ( S / m ) * rhoA * VA * ( x2 / ( r * H ) ) * ( x4 + thetadot * x2 ) - ...
          (1/2) * Cd * ( S / m ) * ( rhoA / VA ) * ( thetadot * x4 + thetadot^2 * x2 ) * ( x4 + thetadot * x2 ) - ...
          (1/2) * Cd * ( S / m ) * rhoA * VA * thetadot;

A(4,3) = 3 * ( mu * x1 * x3 / r^5 ) * ( 1 - (5/2) * J2 * ( Re / r )^2 * ( 7 * ( x3 / r )^2 - 3 ) ) + ...
          (1/2) * Cd * ( S / m ) * rhoA * VA * ( x3 / ( r * H ) ) * ( x4 + thetadot * x2 );

A(4,4) = -(1/2) * Cd * ( S / m ) * ( rhoA / VA ) * ( x4 + thetadot * x2 )^2 - ...
          (1/2) * Cd * ( S / m ) * rhoA * VA;

A(4,5) = -(1/2) * Cd * ( S / m ) * ( rhoA / VA ) * ( x5 - thetadot * x1 ) * ( x4 + thetadot * x2 );
A(4,6) = -(1/2) * Cd * ( S / m ) * ( rhoA / VA ) * x6 * ( x4 + thetadot * x2 );
A(4,7) = -( x1 / r^3 ) * ( 1 - (3/2) * J2 * ( Re / r )^2 * ( 5 * ( x3 / r )^2 - 1 ) );
A(4,8) = (3/2) * ( mu * x1 / r^3 ) * ( Re / r )^2 * ( 5 * ( x3 / r )^2 - 1 );
A(4,9) = -(1/2) * ( S / m ) * rhoA * VA * ( x4 + thetadot * x2 );

A(5,1) = 3 * ( mu * x1 * x2 / r^5 ) * ( 1 - (5/2) * J2 * ( Re / r )^2 * ( 7 * ( x3 / r )^2 - 1 ) ) + ...
          (1/2) * Cd * ( S / m ) * rhoA * VA * ( x1 / ( r * H ) ) * ( x5 - thetadot * x1 ) - ...
          (1/2) * Cd * ( S / m ) * ( rhoA / VA ) * ( thetadot^2 * x1 - thetadot * x5 ) * ( x5 - thetadot * x1 ) + ...
          (1/2) * Cd * ( S / m ) * rhoA * VA * thetadot;

A(5,2) = -( mu / r^3 ) * ( 1 - (3/2) * J2 * ( Re / r )^2 * ( 5 * ( x3 / r )^2 - 1 ) ) + ...
          3 * ( mu * x2^2 / r^5 ) * ( 1 - (5/2) * J2 * ( Re / r )^2 * ( 7 * ( x3 / r )^2 - 1 ) ) + ...
          (1/2) * Cd * ( S / m ) * rhoA * VA * ( x2 / ( r * H ) ) * ( x5 - thetadot * x1 ) - ...
          (1/2) * Cd * ( S / m ) * ( rhoA / VA ) * ( thetadot * x4 + thetadot^2 * x2 ) * ( x5 - thetadot * x1 );

A(5,3) = 3 * ( mu * x2 * x3 / r^5 ) * ( 1 - (5/2) * J2 * ( Re / r )^2 * ( 7 * ( x3 / r )^2 - 3 ) ) + ...
          (1/2) * Cd * ( S / m ) * rhoA * VA * ( x3 / ( r * H ) ) * ( x5 - thetadot * x1 );

A(5,4) = -(1/2) * Cd * ( S / m ) * ( rhoA / VA ) * ( x5 - thetadot * x1 ) * ( x4 + thetadot * x2 );
A(5,5) = -(1/2) * Cd * ( S / m ) * ( rhoA / VA ) * ( x5 - thetadot * x1 )^2 - ...
          (1/2) * Cd * ( S / m ) * rhoA * VA;
A(5,6) = -(1/2) * Cd * ( S / m ) * ( rhoA / VA ) * x6 * ( x5 - thetadot * x1 );
A(5,7) = -( x2 / r^3 ) * ( 1 - (3/2) * J2 * ( Re / r )^2 * ( 5 * ( x3 / r )^2 - 1 ) );
A(5,8) = (3/2) * ( mu * x2 / r^3 ) * ( Re / r )^2 * ( 5 * ( x3 / r )^2 - 1 );
A(5,9) = -(1/2) * ( S / m ) * rhoA * VA * ( x5 - thetadot * x1 );

A(6,1) = 3 * ( mu * x1 * x3 / r^5 ) * ( 1 - (5/2) * J2 * ( Re / r )^2 * ( 7 * ( x3 / r )^2 - 3 ) ) + ...
          (1/2) * Cd * ( S / m ) * rhoA * VA * ( x6 * x1 / ( r * H ) ) - ...
          (1/2) * Cd * ( S / m ) * ( rhoA / VA ) * x6 * ( thetadot^2 * x1 - thetadot * x5 );

A(6,2) = 3 * ( mu * x2 * x3 / r^5 ) * ( 1 - (5/2) * J2 * ( Re / r )^2 * ( 7 * ( x3 / r )^2 - 3 ) ) + ...
          (1/2) * Cd * ( S / m ) * rhoA * VA * ( x6 * x2 / ( r * H ) ) - ...
          (1/2) * Cd * ( S / m ) * ( rhoA / VA ) * x6 * ( thetadot * x4 + thetadot^2 * x2 );

A(6,3) = -( mu / r^3 ) * ( 1 - (3/2) * J2 * ( Re / r )^2 * ( 5 * ( x3 / r )^2 - 3 ) ) + ...
          3 * ( mu * x3^2 / r^5 ) * ( 1 - (5/2) * J2 * ( Re / r )^2 * ( 7 * ( x3 / r )^2 - 5 ) ) + ...
          (1/2) * Cd * ( S / m ) * rhoA * VA * ( x3 * x6 / ( r * H ) );

A(6,4) = -(1/2) * Cd * ( S / m ) * ( rhoA / VA ) * x6 * ( x4 + thetadot * x2 );
A(6,5) = -(1/2) * Cd * ( S / m ) * ( rhoA / VA ) * x6 * ( x5 - thetadot * x1 );
A(6,6) = -(1/2) * Cd * ( S / m ) * ( rhoA / VA ) * x6^2 - ...
          (1/2) * Cd * ( S / m ) * rhoA * VA;
A(6,7) = -( x3 / r^3 ) * ( 1 - (3/2) * J2 * ( Re / r )^2 * ( 5 * ( x3 / r )^2 - 3 ) );
A(6,8) = (3/2) * ( mu * x3 / r^3 ) * ( Re / r )^2 * ( 5 * ( x3 / r )^2 - 3 );
A(6,9) = -(1/2) * ( S / m ) * rhoA * VA * x6;

% Rows 7:9 remain zero because mu, J2, and CD are modeled as constants.

end
