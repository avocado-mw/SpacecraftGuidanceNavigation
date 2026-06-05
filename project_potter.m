%-------------------------------------------------------------------------%
%    project_potter
%
%  DESCRIPTION:
%     Potter square-root sequential filter for the MAE 182 project.
%
%  MAE 182 -- Spacecraft Guidance and Navigation
%-------------------------------------------------------------------------%
clearvars; close all; clc

%% 0) Shared setup
S = project_setup();

%% 1) Initial conditions
x_hat_prev = S.x0_bar;
W_prev     = chol( S.P0_bar , 'lower' );
Xref0_ode  = S.X0_ode;
Phi0       = S.Phi0;

nObs = size( S.Y , 1 );

x_hat    = zeros( S.nState , nObs );
P_store  = zeros( S.nState , S.nState , nObs );
pre_fit  = zeros( 2 , nObs );
post_fit = zeros( 2 , nObs );

t_prev = 0;
t_obs  = S.Yraw(:,1)';

for k = 1:nObs

    t_k = S.Yraw(k,1);

    [ x_ref , Phi_k ] = propagate_step( S , Xref0_ode , t_prev , t_k );

    x_bar = Phi_k * x_hat_prev;
    W_bar = Phi_k * W_prev;

    [ Htilde , rho , rhodot ] = computed_Htilde( t_k , x_ref , S.Y(k,2) );
    y_k = S.Y(k,3:4)' - [ rho ; rhodot ];

    pre_fit(:,k) = y_k - Htilde * x_bar;
    x_upd = x_bar;
    W_upd = W_bar;

    innov1 = y_k(1) - Htilde(1,:) * x_upd;
    [ W_upd , x_upd ] = potter_scalar_update( W_upd , x_upd , Htilde(1,:) , innov1 , S.R(1,1) );

    innov2 = y_k(2) - Htilde(2,:) * x_upd;
    [ W_upd , x_upd ] = potter_scalar_update( W_upd , x_upd , Htilde(2,:) , innov2 , S.R(2,2) );

    x_hat(:,k) = x_upd;
    P_k = W_upd * W_upd';

    if S.symmetrizeCov
        P_k = 0.5 * ( P_k + P_k' );
    end

    post_fit(:,k) = y_k - Htilde * x_hat(:,k);
    P_store(:,:,k) = P_k;

    x_hat_prev = x_hat(:,k);
    W_prev     = W_upd;
    Xref0_ode  = [ x_ref ; Phi0 ];
    t_prev     = t_k;

end

%% 2) RMS values and plots
rho_rms_pre    = sqrt( mean( pre_fit(1,:).^2 ) );
rho_rms_post   = sqrt( mean( post_fit(1,:).^2 ) );
rhodot_rms_pre = sqrt( mean( pre_fit(2,:).^2 ) );
rhodot_rms_post = sqrt( mean( post_fit(2,:).^2 ) );

fprintf( '\nPotter square-root filter summary:\n' );
fprintf( '  Pre-fit  range RMS      = %.4f m\n'   , rho_rms_pre );
fprintf( '  Post-fit range RMS     = %.4f m\n'   , rho_rms_post );
fprintf( '  Pre-fit  range-rate RMS = %.6f m/s\n', rhodot_rms_pre );
fprintf( '  Post-fit range-rate RMS = %.6f m/s\n', rhodot_rms_post );

if S.makePlots
    plot_filter_residuals( t_obs , pre_fit , post_fit , 'Potter' );

    figure( 'Name' , 'Potter Position Error Ellipsoid' );
    Ppos = P_store(1:3,1:3,end);
    [ Rell , semi ] = position_ellipsoid_semiaxes( Ppos );
    plotEllipsoid( Rell , 3 * semi );
    title( 'Potter Position Error Ellipsoid (3-sigma, final epoch)' );
    xlabel( 'x [m]' ); ylabel( 'y [m]' ); zlabel( 'z [m]' );
end

save project_potter_results x_hat P_store pre_fit post_fit ...
     rho_rms_pre rho_rms_post rhodot_rms_pre rhodot_rms_post
