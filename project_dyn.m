%-------------------------------------------------------------------------%
%    project_dyn
%
%  DESCRIPTION:
%     Dynamics model used in MAE 182 Project -- Part I.
%
%     The propagated dynamic subset is
%         [ r ; v ; mu ; J2 ; CD ]
%     together with the 9x9 state transition matrix Phi.
%
%     The acceleration model includes:
%       - two-body gravity,
%       - J2 oblateness,
%       - atmospheric drag with an exponential density model.
%
%     The final twelve stored state components (station coordinates) are
%     treated as constants during this propagation and therefore have zero
%     dynamics here.
%
%  INPUT:
%     x = [ r; v; mu; J2; CD; station coordinates; reshape(Phi,81,1) ]
%
%  OUTPUT:
%     dx = time derivative of the full propagated state
%
%  MAE 182 -- Spacecraft Guidance and Navigation
%-------------------------------------------------------------------------%
function dx = project_dyn( ~ , x )

%% 1) Constants and estimated parameters
Re = 6378136.3;         % Earth radius [m]
mu = x(7);              % Earth gravitational parameter [m^3/s^2]
J2 = x(8);              % Earth oblateness coefficient [-]
Cd = x(9);              % Drag coefficient [-]

S        = 3.0;         % spacecraft projected area [m^2]
m        = 970;         % spacecraft mass [kg]
rho0     = 3.614e-13;   % reference atmospheric density [kg/m^3]
r0       = 7e5 + Re;    % reference radius for density model [m]
H        = 88667.0;     % scale height [m]
thetadot = 7.2921158553e-5;  % Earth rotation rate [rad/s]

%% 2) Satellite position and velocity
rvec = x(1:3);
vvec = x(4:6);
r    = norm( rvec );

%% 3) Two-body gravity
a_g_2B = -mu * rvec / r^3;

%% 4) J2 perturbation
z_over_r = x(3) / r;
commonJ2 = ( 3 / 2 ) * J2 * ( Re / r )^2;

a_g_J2      = zeros( 3 , 1 );
a_g_J2(1)   = ( mu * x(1) / r^3 ) * commonJ2 * ( 5 * z_over_r^2 - 1 );
a_g_J2(2)   = ( mu * x(2) / r^3 ) * commonJ2 * ( 5 * z_over_r^2 - 1 );
a_g_J2(3)   = ( mu * x(3) / r^3 ) * commonJ2 * ( 5 * z_over_r^2 - 3 );

%% 5) Atmospheric drag
rhoA = rho0 * exp( -( r - r0 ) / H );
Va   = [ x(4) + thetadot * x(2) ; ...
         x(5) - thetadot * x(1) ; ...
         x(6) ];

a_drag = -( 1 / 2 ) * Cd * ( S / m ) * rhoA * norm( Va ) * Va;

%% 6) Total translational acceleration
a_tot = a_g_2B + a_g_J2 + a_drag;

%% 7) STM propagation
Phi    = reshape( x(19:end) , 9 , 9 );
A      = computed_A( x );
Phidot = A * Phi;
Phidot = reshape( Phidot , 9 * 9 , 1 );

%% 8) Assemble the full state derivative
% Station coordinates are constants in this Part I setup.
dx = [ vvec ; a_tot ; zeros( 12 , 1 ) ; Phidot ];

end
