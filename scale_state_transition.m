%-------------------------------------------------------------------------%
%    scale_state_transition
%
%  DESCRIPTION:
%     Maps STM to normalized coordinates with x = S*z:
%         Phi_z = S^{-1} Phi S
%
%  MAE 182 -- Spacecraft Guidance and Navigation
%-------------------------------------------------------------------------%
function Phiz = scale_state_transition( S , Phi )

Sinv = diag( 1 ./ diag( S ) );
Phiz = Sinv * Phi * S;

end
