%-------------------------------------------------------------------------%
%    position_covariance_metrics
%
%  DESCRIPTION:
%     Summarizes the 3x3 position covariance block.
%
%  MAE 182 -- Spacecraft Guidance and Navigation
%-------------------------------------------------------------------------%
function m = position_covariance_metrics( P )

Ppos = real( 0.5 * ( P(1:3,1:3) + P(1:3,1:3)' ) );
[ ~ , D ] = eig( Ppos );
eigvals = max( real( diag( D ) ) , 0 );

m.sigma_xyz   = sqrt( max( real( diag( Ppos ) ) , 0 ) );
m.sigma_3axis = sqrt( eigvals );
m.sigma_rms   = sqrt( mean( m.sigma_3axis.^2 ) );
m.trace       = trace( Ppos );
m.det         = det( Ppos );

end
