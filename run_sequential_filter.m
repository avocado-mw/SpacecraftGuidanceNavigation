%-------------------------------------------------------------------------%
%    run_sequential_filter
%
%  DESCRIPTION:
%     Runs the Joseph-form sequential filter in normalized state coordinates
%     (x = S*z) for numerical conditioning.  Returns summary metrics.
%
%  MAE 182 -- Spacecraft Guidance and Navigation
%-------------------------------------------------------------------------%
function out = run_sequential_filter( S , P0_bar )

nObs    = size( S.Y , 1 );
nState  = S.nState;
I18     = eye( nState );
S_scale = build_state_scaling();

x_hat_prev = S_scale \ S.x0_bar;
P_prev     = scale_state_covariance( S_scale , P0_bar );
Xref0_ode  = S.X0_ode;
Phi0       = S.Phi0;
t_prev     = 0;

pre_fit  = zeros( 2 , nObs );
post_fit = zeros( 2 , nObs );

for k = 1:nObs

    t_k = S.Yraw(k,1);

    [ x_ref , Phi_k ] = propagate_step( S , Xref0_ode , t_prev , t_k );
    Phi_k = scale_state_transition( S_scale , Phi_k );

    x_bar = Phi_k * x_hat_prev;
    P_bar = Phi_k * P_prev * Phi_k';

    if S.symmetrizeCov
        P_bar = 0.5 * ( P_bar + P_bar' );
    end

    [ Htilde , rho , rhodot ] = computed_Htilde( t_k , x_ref , S.Y(k,2) );
    H_n = Htilde * S_scale;
    y_k = S.Y(k,3:4)' - [ rho ; rhodot ];

    Syy = H_n * P_bar * H_n' + S.R;
    K   = P_bar * H_n' / Syy;

    pre_fit(:,k)  = y_k - Htilde * ( S_scale * x_bar );
    x_hat_k       = x_bar + K * ( y_k - H_n * x_bar );
    P_k           = sequential_covariance_update( S , P_bar , K , H_n , I18 );
    post_fit(:,k) = y_k - Htilde * ( S_scale * x_hat_k );

    x_hat_prev = x_hat_k;
    P_prev     = P_k;
    Xref0_ode  = [ x_ref ; Phi0 ];
    t_prev     = t_k;

end

P_final = unscale_state_covariance( S_scale , P_k );

out.P_final         = P_final;
out.rho_rms_pre     = sqrt( mean( pre_fit(1,:).^2 ) );
out.rhodot_rms_pre  = sqrt( mean( pre_fit(2,:).^2 ) );
out.rho_rms_post    = sqrt( mean( post_fit(1,:).^2 ) );
out.rhodot_rms_post = sqrt( mean( post_fit(2,:).^2 ) );
out.pos_metrics     = position_covariance_metrics( P_final );

end
