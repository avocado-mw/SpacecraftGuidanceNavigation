%-------------------------------------------------------------------------%
%    run_batch_filter
%
%  DESCRIPTION:
%     Runs the batch least-squares processor with state normalization and
%     QR solution of the normal equations.  Returns summary metrics.
%
%  MAE 182 -- Spacecraft Guidance and Navigation
%-------------------------------------------------------------------------%
function out = run_batch_filter( S , P0_bar )

nState  = S.nState;
nObs    = size( S.Y , 1 );
nIter   = S.batch_nIter;
S_scale = build_state_scaling();

X0_ode = S.X0_ode;
x0_bar = S.x0_bar;

for pass = 1:nIter

    [ t20 , x20 , Phi20 ] = propagate_full_arc( S , X0_ode );

    H_stack = zeros( 2 * nObs , nState );
    y_stack = zeros( 2 * nObs , 1 );

    for k = 1:nObs
        j = S.Yidx(k,1);
        rows = ( 2 * ( k - 1 ) + 1 ) : ( 2 * k );

        [ Htilde , rho , rhodot ] = computed_Htilde( t20(j) , x20(:,j) , S.Y(k,2) );
        H_stack(rows,:) = Htilde * Phi20(:,:,j);
        y_stack(rows)   = S.Y(k,3:4)' - [ rho ; rhodot ];
    end

    y = reshape( y_stack , 2 , nObs );

    [ x0_hat , P0 ] = batch_solve_qr( S , P0_bar , H_stack , y_stack , x0_bar , S_scale );

    X0_ode(1:nState) = X0_ode(1:nState) + x0_hat;
    x0_bar           = zeros( nState , 1 );

    if pass == 1
        rho_rms_pre    = sqrt( mean( y(1,:).^2 ) );
        rhodot_rms_pre = sqrt( mean( y(2,:).^2 ) );
    end

    rho_rms_post    = sqrt( mean( y(1,:).^2 ) );
    rhodot_rms_post = sqrt( mean( y(2,:).^2 ) );

end

out.P_final         = P0;
out.rho_rms_pre     = rho_rms_pre;
out.rhodot_rms_pre  = rhodot_rms_pre;
out.rho_rms_post    = rho_rms_post;
out.rhodot_rms_post = rhodot_rms_post;
out.pos_metrics     = position_covariance_metrics( P0 );

end
