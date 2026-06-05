%-------------------------------------------------------------------------%
%    scale_state_covariance
%
%  DESCRIPTION:
%     Maps physical covariance P to normalized coordinates with x = S*z:
%         P_z = S^{-1} P S^{-1}
%
%  MAE 182 -- Spacecraft Guidance and Navigation
%-------------------------------------------------------------------------%
function Pz = scale_state_covariance( S , P )

Sinv = diag( 1 ./ diag( S ) );
Pz   = Sinv * P * Sinv;

end
