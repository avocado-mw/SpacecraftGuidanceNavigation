%-------------------------------------------------------------------------%
%    project_part_1
%
%  DESCRIPTION:
%     Cleaned reference solution for MAE 182 Project -- Part I.
%
%     This script follows the same overall organization as the learner
%     template, but with the required computations completed.  The workflow is:
%       (1) define the a priori state and STM initial conditions,
%       (2) propagate the reference trajectory and the 9x9 STM,
%       (3) resample the solution every 20 seconds,
%       (4) compute predicted range / range-rate and Htilde,
%       (5) form residuals and RMS values,
%       (6) report the requested Part I outputs.
%
%  FILES EXPECTED IN THE SAME FOLDER:
%     project_dyn.m
%     computed_A.m
%     computed_Htilde.m
%     computed_obs.m
%     project_obs_data.txt
%
%  NOTES:
%     - The propagated dynamic state is the 9-state subset
%           [ r ; v ; mu ; J2 ; CD ]
%       while the station coordinates are appended to the stored state but are
%       treated as constants during this Part I propagation.
%     - The STM propagated by project_dyn.m is therefore 9x9.  When storing the
%       full 18x18 Phi matrix required by the linearized observation model, the
%       lower-right station-coordinate block is augmented with an identity.
%
%  MAE 182 -- Spacecraft Guidance and Navigation
%-------------------------------------------------------------------------%
clearvars; close all; clc

%% 0) A priori / initial values from the handout
Phi0 = reshape( eye( 9 ) , 9 * 9 , 1 );

% Approximate initial satellite state in ECI [m; m/s]
xv0 = [ [ 757700  , 5222607 , 4851500 ]' ; ...
        [ 2213.21 , 4678.34 , -5371.30 ]' ];

% Tracking-station coordinates in ECEF [m]
Xsite0 = [ [ -5127510 , -3794160 ,       0 ]' ; ...
           [  3860910 ,  3238490 , 3898094 ]' ; ...
           [   549505 , -1380872 , 6182197 ]' ];

% Approximate physical parameters [mu; J2; CD]
param0 = [ 3.986004415e14 ; 1.082626925638815e-3 ; 2.0 ];

% Stored state vector used by the observation model:
%   x = [r; v; mu; J2; CD; station coordinates]
X0_state = [ xv0 ; param0 ; Xsite0 ];

% Full ODE state includes the reshaped 9x9 STM appended at the end.
X0 = [ X0_state ; Phi0 ];

%% 1) Integrate the reference trajectory and the 9x9 STM
% The assignment specifies propagation from t = 0 to t = 18340 s.
tol     = 1e-12;
tstep   = 0.1;
time    = 0:tstep:18340;
options = odeset( 'RelTol' , tol , 'AbsTol' , tol );

[ ~ , X ] = ode45( @project_dyn , time , X0 , options );

%% 2) Unpack the propagated state and the 9x9 STM history
% First 18 entries are [r; v; mu; J2; CD; station coordinates].
x_temp = X(:,1:18)';

% Quick sanity check: the reference trajectory should resemble a near-polar
% low-Earth orbit.
figure(1)
plot3( x_temp(1,:) , x_temp(2,:) , x_temp(3,:) , 'k-' )
grid on; axis equal
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]')
title('Reference trajectory in ECI')

% Reshape the propagated STM blocks.
Phi_temp = X(:,19:end)';
Phi_temp = reshape( Phi_temp , 9 , 9 , length(time) );

%% 3) Resample the propagated solution every 20 seconds
% Observation data are provided every 20 seconds, so we store a matching grid.
rsmpl     = round( 20 / tstep );
lng_rsmpl = round( length(time) / rsmpl ) + 1;

t   = zeros( lng_rsmpl , 1 );
x   = zeros( 18 , lng_rsmpl );
Phi = zeros( 18 , 18 , lng_rsmpl );

for k = 1:lng_rsmpl
    idx    = rsmpl * ( k - 1 ) + 1;

    t(k)   = time(idx);
    x(:,k) = x_temp(:,idx);

    % The propagated STM covers only the 9 dynamic states [r; v; mu; J2; CD].
    % Station coordinates are constant in this Part I reference propagation, so
    % the corresponding STM block is the identity.
    Phi(:,:,k) = [ Phi_temp(:,:,idx) , zeros( 9 , 9 ) ; ...
                   zeros( 9 , 9 )    , eye( 9 , 9 )   ];
end

%% 4) Load and organize the observation data
Y = load( 'project_obs_data.txt' );

% Convert time [s] to an index into the 20-second resampled arrays.
Y(:,1) = round( Y(:,1) / 20 ) + 1;

% Map station IDs {101,337,394} to local indices {1,2,3} used by the
% helper functions.
Y( Y(:,2) == 101 , 2 ) = 1;
Y( Y(:,2) == 337 , 2 ) = 2;
Y( Y(:,2) == 394 , 2 ) = 3;

lng_obs = size( Y , 1 );

%% 5) Compute predicted observations, Htilde, and residuals
rho    = zeros( lng_obs , 1 );
rhodot = zeros( lng_obs , 1 );
Htilde = zeros( 2 , 18 , lng_obs );
y      = zeros( 2 , lng_obs );

for k = 1:lng_obs
    j = Y(k,1);    % resampled time index

    [ Htilde(:,:,k) , rho(k) , rhodot(k) ] = computed_Htilde( t(j) , x(:,j) , Y(k,2) );

    % Residuals are stored in the standard O - C form.
    y(:,k) = Y(k,3:4)' - [ rho(k) ; rhodot(k) ];
end

%% 6) Compute RMS values for the residuals
rho_res    = y(1,:);
rhodot_res = y(2,:);

rho_rms    = sqrt( sum( rho_res.^2 )    / lng_obs );
rhodot_rms = sqrt( sum( rhodot_res.^2 ) / lng_obs );

disp( [ 'range_rms = '      , num2str( rho_rms )    , ' m'   ] );
disp( [ 'range_rate_rms = ' , num2str( rhodot_rms ) , ' m/s' ] );

%% 7) Residual plots
figure(2)
plot( rho_res , 'k-' )
grid on
hTitle  = title( 'Range Residuals' );
hXLabel = xlabel( 'observation number' );
hYLabel = ylabel( 'O - C [m]' );
set( gca , 'FontName' , 'Helvetica' );
set( [ hTitle , hXLabel , hYLabel ] , 'FontName' , 'AvantGarde' );
set( [ hTitle , hXLabel , hYLabel ] , 'FontSize' , 12 );

figure(3)
plot( rhodot_res , 'k-' )
grid on
hTitle  = title( 'Range-Rate Residuals' );
hXLabel = xlabel( 'observation number' );
hYLabel = ylabel( 'O - C [m/s]' );
set( gca , 'FontName' , 'Helvetica' );
set( [ hTitle , hXLabel , hYLabel ] , 'FontName' , 'AvantGarde' );
set( [ hTitle , hXLabel , hYLabel ] , 'FontSize' , 12 );

%% 8) Requested outputs for Part I
% The assignment asks for Htilde(t0), the upper-left 6x9 block of A(t0), and
% the propagated r(tf), v(tf), and Phi(tf,t0).
A0         = computed_A( X0_state );
Htilde0    = Htilde(:,:,1);
xvf        = x(1:6,end);
Phif       = Phi(:,:,end);
residuals0 = y(:,1);

save project_part_1_results A0 Htilde0 xvf Phif residuals0 rho_rms rhodot_rms

%% 9) Reminder of what to report
fprintf('\nProject Part I outputs:\n');
fprintf('  • Htilde(t0)\n');
fprintf('  • A(t0) upper-left 6x9 block\n');
fprintf('  • r(tf), v(tf), and Phi(tf,t0)\n');
fprintf('  • range and range-rate RMS values\n');
