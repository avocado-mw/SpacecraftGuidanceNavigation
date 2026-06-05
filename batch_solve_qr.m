%-------------------------------------------------------------------------%
%    batch_solve_qr
%
%  DESCRIPTION:
%     Solves the batch normal equations using state normalization and QR
%     decomposition instead of explicit matrix inversion.
%
%     Minimizes  ||R^{-1/2}(y - Hx)||^2 + ||P0^{-1/2}(x - x0)||^2
%     with x = S*z and augmented system  A z ≈ b  solved via QR.
%
%  MAE 182 -- Spacecraft Guidance and Navigation
%-------------------------------------------------------------------------%
function [ x_hat , P_post ] = batch_solve_qr( S , P0_bar , H_stack , y_stack , x0_bar , S_scale )

nState = size( H_stack , 2 );
nObs   = size( y_stack , 1 ) / 2;

Rw = diag( [ 1 / sqrt( S.R(1,1) ) , 1 / sqrt( S.R(2,2) ) ] );

H_w = zeros( size( H_stack ) );
y_w = zeros( size( y_stack ) );

for k = 1:nObs
    rows      = ( 2 * ( k - 1 ) + 1 ) : ( 2 * k );
    H_w(rows,:) = Rw * H_stack(rows,:);
    y_w(rows)   = Rw * y_stack(rows);
end

H_n   = H_w * S_scale;
P0_n  = scale_state_covariance( S_scale , P0_bar );
z0    = S_scale \ x0_bar;

prior_sqrt = diag( 1 ./ sqrt( max( diag( P0_n ) , eps ) ) );
A          = [ H_n ; prior_sqrt ];
b          = [ y_w ; prior_sqrt * z0 ];

[ Q , R_qr ] = qr( A , 0 );
z_hat        = R_qr \ ( Q' * b );
x_hat        = S_scale * z_hat;

P_z    = ( R_qr' * R_qr ) \ eye( nState );
P_post = unscale_state_covariance( S_scale , P_z );

if S.symmetrizeCov
    P_post = 0.5 * ( P_post + P_post' );
end

P_post = real( P_post );

end
