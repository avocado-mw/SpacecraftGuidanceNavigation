%-------------------------------------------------------------------------%
%    project_batch
%
%  DESCRIPTION:
%     Batch least-squares processor for the MAE 182 orbit-determination
%     project.
%
%  MAE 182 -- Spacecraft Guidance and Navigation
%-------------------------------------------------------------------------%
clearvars; close all; clc

%% 0) Shared setup
S = project_setup();
nState = S.nState;

%% 1) Initialize the batch iteration
X0_ode = S.X0_ode;
x0_bar = S.x0_bar;
P0_bar = S.P0_bar;

nIter           = S.batch_nIter;
x0_hat_hist     = zeros( nState , nIter );
rho_rms_hist    = zeros( nIter , 1 );
rhodot_rms_hist = zeros( nIter , 1 );

nObs = size( S.Y , 1 );
y_all = zeros( 2 , nObs , nIter );

for pass = 1:nIter

    %% 2) Propagate the reference trajectory and resampled STM history
    [ t20 , x20 , Phi20 ] = propagate_full_arc( S , X0_ode );

    %% 3) Stack the measurement model for QR solve
    S_scale = build_state_scaling();
    H_stack = zeros( 2 * nObs , nState );
    y_stack = zeros( 2 * nObs , 1 );

    rho    = zeros( nObs , 1 );
    rhodot = zeros( nObs , 1 );
    y      = zeros( 2 , nObs );

    %% 4) Loop over the observations
    for k = 1:nObs
        j = S.Yidx(k,1);
        rows = ( 2 * ( k - 1 ) + 1 ) : ( 2 * k );

        [ Htilde , rho(k) , rhodot(k) ] = computed_Htilde( t20(j) , x20(:,j) , S.Y(k,2) );

        H_stack(rows,:) = Htilde * Phi20(:,:,j);
        y(:,k)          = S.Y(k,3:4)' - [ rho(k) ; rhodot(k) ];
        y_stack(rows)   = y(:,k);
    end

    y_all(:,:,pass) = y;

    %% 5) RMS values for this pass
    rho_rms_hist(pass)    = sqrt( mean( y(1,:).^2 ) );
    rhodot_rms_hist(pass) = sqrt( mean( y(2,:).^2 ) );

    disp( [ 'Pass ' , num2str(pass) , ': range_rms = ' , num2str( rho_rms_hist(pass) ) , ...
            ' m, range_rate_rms = ' , num2str( rhodot_rms_hist(pass) ) , ' m/s' ] );

    %% 6) QR solve in normalized coordinates and update the nominal state
    [ x0_hat , P0 ] = batch_solve_qr( S , P0_bar , H_stack , y_stack , x0_bar , S_scale );

    X0_ode(1:nState) = X0_ode(1:nState) + x0_hat;
    x0_hat_hist(:,pass) = x0_hat;

    % Deviation from the updated nominal is zero at the start of the next pass.
    x0_bar = zeros( nState , 1 );

end

%% 7) Summary outputs and plots
pre_fit  = y_all(:,:,1);
post_fit = y_all(:,:,nIter);
P_final  = P0;

rho_rms_pre     = rho_rms_hist(1);
rho_rms_post    = rho_rms_hist(end);
rhodot_rms_pre  = rhodot_rms_hist(1);
rhodot_rms_post = rhodot_rms_hist(end);

fprintf( '\nBatch least squares summary:\n' );
fprintf( '  Measurement noise R = diag([%.2e, %.2e]) [m^2, (m/s)^2]\n' , S.R(1,1) , S.R(2,2) );
fprintf( '  Pre-fit  range RMS      = %.4f m\n'   , rho_rms_pre );
fprintf( '  Post-fit range RMS     = %.4f m\n'   , rho_rms_post );
fprintf( '  Pre-fit  range-rate RMS = %.6f m/s\n', rhodot_rms_pre );
fprintf( '  Post-fit range-rate RMS = %.6f m/s\n', rhodot_rms_post );

pos_sigma = sqrt( max( real( diag( P_final(1:3,1:3) ) ) , 0 ) );
fprintf( '  Position 1-sigma (x,y,z) = [%.3f, %.3f, %.3f] m\n' , pos_sigma );

if S.makePlots
    t_obs = S.Yraw(:,1)';
    plot_filter_residuals( t_obs , pre_fit , post_fit , 'Batch' );

    figure( 'Name' , 'Batch Position Error Ellipsoid' );
    Ppos = P_final(1:3,1:3);
    [ Rell , semi ] = position_ellipsoid_semiaxes( Ppos );
    plotEllipsoid( Rell , 3 * semi );
    title( 'Batch Position Error Ellipsoid (3-sigma)' );
    xlabel( 'x [m]' ); ylabel( 'y [m]' ); zlabel( 'z [m]' );
end

save project_batch_results x0_hat_hist rho_rms_hist rhodot_rms_hist ...
     pre_fit post_fit P_final rho_rms_pre rho_rms_post rhodot_rms_pre rhodot_rms_post
