%-------------------------------------------------------------------------%
%    run_batch_filter
%
%  DESCRIPTION:
%     Runs the batch least-squares processor and returns summary metrics.
%
%  MAE 182 -- Spacecraft Guidance and Navigation
%-------------------------------------------------------------------------%
function out = run_batch_filter( S , P0_bar )

nState = S.nState;
nObs   = size( S.Y , 1 );
nIter  = S.batch_nIter;

X0_ode = S.X0_ode;
x0_bar = S.x0_bar;

for pass = 1:nIter

    [ t20 , x20 , Phi20 ] = propagate_full_arc( S , X0_ode );

    Lambda = P0_bar \ eye( nState );
    N      = Lambda * x0_bar;

    y = zeros( 2 , nObs );

    for k = 1:nObs
        j = S.Yidx(k,1);
        [ Htilde , rho , rhodot ] = computed_Htilde( t20(j) , x20(:,j) , S.Y(k,2) );
        Hk = Htilde * Phi20(:,:,j);
        y(:,k) = S.Y(k,3:4)' - [ rho ; rhodot ];
        Lambda = Lambda + Hk' * ( S.R \ Hk );
        N      = N      + Hk' * ( S.R \ y(:,k) );
    end

    x0_hat = Lambda \ N;
    P0     = inv( Lambda );

    if S.symmetrizeCov
        P0 = 0.5 * ( P0 + P0' );
    end

    X0_ode(1:nState) = X0_ode(1:nState) + x0_hat;
    x0_bar = zeros( nState , 1 );

    if pass == 1
        rho_rms_pre    = sqrt( mean( y(1,:).^2 ) );
        rhodot_rms_pre = sqrt( mean( y(2,:).^2 ) );
    end

    rho_rms_post    = sqrt( mean( y(1,:).^2 ) );
    rhodot_rms_post = sqrt( mean( y(2,:).^2 ) );

end

out.P_final          = P0;
out.rho_rms_pre      = rho_rms_pre;
out.rhodot_rms_pre   = rhodot_rms_pre;
out.rho_rms_post     = rho_rms_post;
out.rhodot_rms_post  = rhodot_rms_post;
out.pos_metrics      = position_covariance_metrics( P0 );

end
