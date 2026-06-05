%-------------------------------------------------------------------------%
%    position_ellipsoid_semiaxes
%
%  DESCRIPTION:
%     Returns real, non-negative semi-axis lengths for a 3x3 position
%     covariance block.  Negative eigenvalues from round-off are clipped
%     before taking square roots so plotEllipsoid never receives complex
%     coordinates.
%
%  MAE 182 -- Spacecraft Guidance and Navigation
%-------------------------------------------------------------------------%
function [ Rell , semi ] = position_ellipsoid_semiaxes( Ppos )

Ppos = real( 0.5 * ( Ppos + Ppos' ) );
[ Rell , D ] = eig( Ppos );
semi = sqrt( max( real( diag( D ) ) , 0 ) );

end
