%-------------------------------------------------------------------------%
%    project_setup
%
%  DESCRIPTION:
%     Shared configuration, constants, initial conditions, and observation
%     data for the MAE 182 orbit-determination project.
%
%  OUTPUT:
%     S  - structure containing constants, observations, and initial values
%
%  MAE 182 -- Spacecraft Guidance and Navigation
%-------------------------------------------------------------------------%
function S = project_setup()

%% 1) ODE and observation settings
S.tol       = 1e-12;
S.tstep     = 0.1;
S.tFinal    = 18340;
S.timeFull  = 0:S.tstep:S.tFinal;
S.rsmpl     = round( 20 / S.tstep );
S.thetadot  = 7.2921158553e-5;

%% 2) Physical constants and approximate parameters
S.Re        = 6378136.3;                % [m]
S.mu0       = 3.986004415e14;           % [m^3/s^2]
S.J20       = 1.082626925638815e-3;     % [-]
S.CD0       = 2.0;                      % [-]
S.area      = 3.0;                      % [m^2]
S.mass      = 970;                      % [kg]
S.rho0      = 3.614e-13;                % [kg/m^3]
S.r0        = 700000.0 + S.Re;          % [m]
S.Hscale    = 88667.0;                  % [m]

%% 3) Approximate initial conditions
S.rv0 = [ [ 757700  , 5222607 , 4851500 ]' ; ...
          [ 2213.21 , 4678.34 , -5371.30 ]' ];

S.Xsite0 = [ [ -5127510 , -3794160 ,       0 ]' ; ...
             [  3860910 ,  3238490 , 3898094 ]' ; ...
             [   549505 , -1380872 , 6182197 ]' ];

S.param0    = [ S.mu0 ; S.J20 ; S.CD0 ];
S.X0_state  = [ S.rv0 ; S.param0 ; S.Xsite0 ];

S.Phi0      = reshape( eye( 9 ) , 9 * 9 , 1 );
S.X0_ode    = [ S.X0_state ; S.Phi0 ];

%% 4) A priori state and covariance for the 18-state estimator
% State ordering:
%   [ r(3); v(3); mu; J2; CD; x_s1(3); x_s2(3); x_s3(3) ]
% Station 101 (s1) is treated as known.  Stations 337 (s2) and 394 (s3) are
% unknown and initialized with sigma = 1 km as specified in the handout.
S.nState   = 18;
S.x0_bar   = zeros( S.nState , 1 );
S.P0_bar   = diag([ ...
    1e6 , 1e6 , 1e6 , ...
    1e6 , 1e6 , 1e6 , ...
    1e20, ...
    1e6 , ...
    1e6 , ...
    1e-10, 1e-10, 1e-10, ...
    1e6 , 1e6 , 1e6 , ...
    1e6 , 1e6 , 1e6 ]);

%% 5) Measurement-noise model  [m] and [m/s]
S.sig_rho    = 0.01;
S.sig_rhodot = 0.001;
S.R          = diag( [ S.sig_rho^2 , S.sig_rhodot^2 ] );
S.W          = inv( S.R );

%% 6) Observation file (all three stations, sorted by time)
S.obsFile = 'project_obs_data.txt';
S.Yraw    = load( S.obsFile );
S.Yraw    = sortrows( S.Yraw , 1 );
S.Y       = S.Yraw;
S.Yidx    = S.Yraw;

S.Y(:,2)    = map_station_ids( S.Y(:,2) );
S.Yidx(:,2) = S.Y(:,2);
S.Yidx(:,1) = round( S.Yraw(:,1) / 20 ) + 1;

%% 7) Default algorithm settings
S.batch_nIter    = 5;
S.saveFigures    = false;
S.makePlots      = true;
S.ekfWarmup      = 100;
S.useJoseph      = true;
S.symmetrizeCov  = true;

end

function stationLocal = map_station_ids( stationIDs )
stationLocal = stationIDs;
stationLocal( stationIDs == 101 ) = 1;
stationLocal( stationIDs == 337 ) = 2;
stationLocal( stationIDs == 394 ) = 3;
end
