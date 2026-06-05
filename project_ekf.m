%-------------------------------------------------------------------------%
%    project_ekf
%
%  DESCRIPTION:
%     Extended Kalman filter for the MAE 182 project.  The reference
%     trajectory is updated after the measurement update once the filter has
%     passed the warmup interval.
%
%  MAE 182 -- Spacecraft Guidance and Navigation
%-------------------------------------------------------------------------%
clearvars; close all; clc

%% 0) Shared setup
S = project_setup( 'sequential' );

fprintf( 'EKF using station ID(s): %s\n' , num2str( S.sequentialStationIDs ) );

%% 1) Initial conditions
x_hat_prev = S.x0_bar;
P_prev     = S.P0_bar;
Xref0_ode  = S.X0_ode;
Phi0       = S.Phi0;
I18        = eye( 18 );

nObs = size( S.Y , 1 );

x_hat    = zeros( 18 , nObs );
P_store  = zeros( 18 , 18 , nObs );
pre_fit  = zeros( 2 , nObs );
post_fit = zeros( 2 , nObs );

t_prev = 0;
t_obs  = S.Yraw(:,1)';

for k = 1:nObs

    t_k = S.Yraw(k,1);
    dt  = t_k - t_prev;

    [ x_hat_prev , P_prev ] = sequential_gap_reset( S , dt , x_hat_prev , P_prev );

    [ x_ref , Phi_k ] = propagate_step( S , Xref0_ode , t_prev , t_k );

    x_bar = Phi_k * x_hat_prev;
    P_bar = Phi_k * P_prev * Phi_k';

    if S.symmetrizeCov
        P_bar = 0.5 * ( P_bar + P_bar' );
    end

    [ Htilde , rho , rhodot ] = computed_Htilde( t_k , x_ref , S.Y(k,2) );
    y_k = S.Y(k,3:4)' - [ rho ; rhodot ];

    if S.useJoseph
        Syy = Htilde * P_bar * Htilde' + S.R;
        K   = P_bar * Htilde' / Syy;
        pre_fit(:,k) = y_k - Htilde * x_bar;
        x_hat(:,k)   = x_bar + K * pre_fit(:,k);
        IKH = I18 - K * Htilde;
        P_k = IKH * P_bar * IKH' + K * S.R * K';
    else
        Syy = Htilde * P_bar * Htilde' + S.R;
        K   = P_bar * Htilde' / Syy;
        pre_fit(:,k) = y_k - Htilde * x_bar;
        x_hat(:,k)   = x_bar + K * pre_fit(:,k);
        P_k = ( I18 - K * Htilde ) * P_bar;
    end

    if S.symmetrizeCov
        P_k = 0.5 * ( P_k + P_k' );
    end

    post_fit(:,k) = y_k - Htilde * x_hat(:,k);
    P_store(:,:,k) = P_k;

    x_ref_updated = x_ref + x_hat(:,k);
    if k >= S.ekfWarmup
        x_hat_prev = zeros( 18 , 1 );
        Xref0_ode  = [ x_ref_updated ; Phi0 ];
    else
        x_hat_prev = x_hat(:,k);
        Xref0_ode  = [ x_ref ; Phi0 ];
    end

    P_prev = P_k;
    t_prev = t_k;

end

%% 2) RMS values and plots
rho_rms_pre    = sqrt( mean( pre_fit(1,:).^2 ) );
rho_rms_post   = sqrt( mean( post_fit(1,:).^2 ) );
rhodot_rms_pre = sqrt( mean( pre_fit(2,:).^2 ) );
rhodot_rms_post = sqrt( mean( post_fit(2,:).^2 ) );

fprintf( '\nExtended Kalman filter summary:\n' );
fprintf( '  Reference updates begin after observation %d\n' , S.ekfWarmup );
fprintf( '  Pre-fit  range RMS      = %.4f m\n'   , rho_rms_pre );
fprintf( '  Post-fit range RMS     = %.4f m\n'   , rho_rms_post );
fprintf( '  Pre-fit  range-rate RMS = %.6f m/s\n', rhodot_rms_pre );
fprintf( '  Post-fit range-rate RMS = %.6f m/s\n', rhodot_rms_post );

if S.makePlots
    plot_filter_residuals( t_obs , pre_fit , post_fit , 'EKF' );

    figure( 'Name' , 'EKF Position Error Ellipsoid' );
    Ppos = P_store(1:3,1:3,end);
    [ Rell , D ] = eig( Ppos );
    semi = sqrt( diag( D ) );
    plotEllipsoid( Rell , semi );
    title( 'EKF Position Error Ellipsoid (final epoch)' );
    xlabel( 'x [m]' ); ylabel( 'y [m]' ); zlabel( 'z [m]' );
end

save project_ekf_results x_hat P_store pre_fit post_fit ...
     rho_rms_pre rho_rms_post rhodot_rms_pre rhodot_rms_post
