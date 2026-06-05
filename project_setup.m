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
function S = project_setup( mode )

if nargin < 1
    mode = 'batch';
end

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

%% 4) A priori state and covariance for the full 18-state estimator
S.x0_bar = zeros( 18 , 1 );
% Station coordinates are known from the handout, so their a priori sigmas are
% set to 1e-5 m (variance 1e-10).  Only the dynamic states are estimated.
S.P0_bar = diag([ ...
    1e6 , 1e6 , 1e6 , ...
    1e6 , 1e6 , 1e6 , ...
    1e20, ...
    1e6 , ...
    1e6 , ...
    1e-10, 1e-10, 1e-10, ...
    1e-10, 1e-10, 1e-10, ...
    1e-10, 1e-10, 1e-10 ]);

%% 5) Measurement-noise model
S.sig_rho    = 0.01;     % [m]
S.sig_rhodot = 0.001;    % [m/s]
S.R          = diag( [ S.sig_rho^2 , S.sig_rhodot^2 ] );
S.W          = diag( [ 1 / S.sig_rho^2 , 1 / S.sig_rhodot^2 ] );

%% 6) Observation file
S.obsFile = 'project_obs_data.txt';
S.YrawAll = load( S.obsFile );

% Batch processing uses all stations.  Sequential filters use a single
% station to avoid mixing incompatible visibility arcs and coordinate
% partials in one filter pass.
S.batchStationIDs      = [ 101 , 337 , 394 ];
S.sequentialStationIDs = 337;

if strcmpi( mode , 'sequential' )
    activeIDs = S.sequentialStationIDs;
else
    activeIDs = S.batchStationIDs;
end

S.Yraw = filter_observations( S.YrawAll , activeIDs );
S.Y    = S.Yraw;
S.Yidx = S.Yraw;

% station IDs --> local indices {1,2,3}
S.Y(:,2)    = map_station_ids( S.Y(:,2) );
S.Yidx(:,2) = S.Y(:,2);

% observation times --> 20-second resampled indices
S.Yidx(:,1) = round( S.Yraw(:,1) / 20 ) + 1;

%% 7) Default algorithm settings
S.batch_nIter    = 3;
S.saveFigures    = false;
S.makePlots      = true;
S.ekfWarmup      = 100;      % update reference after this many observations
S.useJoseph      = true;     % default robust covariance update
S.symmetrizeCov  = true;
S.gapResetTime   = 60;       % reset sequential covariance after [s]

end

function stationLocal = map_station_ids( stationIDs )
stationLocal = stationIDs;
stationLocal( stationIDs == 101 ) = 1;
stationLocal( stationIDs == 337 ) = 2;
stationLocal( stationIDs == 394 ) = 3;
end

function Yfilt = filter_observations( Yall , stationIDs )
keep = ismember( Yall(:,2) , stationIDs );
Yfilt = Yall( keep , : );
Yfilt = sortrows( Yfilt , 1 );
end
