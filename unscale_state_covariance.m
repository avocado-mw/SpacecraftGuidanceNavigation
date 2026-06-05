%-------------------------------------------------------------------------%
%    unscale_state_covariance
%
%  DESCRIPTION:
%     Maps normalized covariance back to physical units with x = S*z:
%         P = S P_z S
%
%  MAE 182 -- Spacecraft Guidance and Navigation
%-------------------------------------------------------------------------%
function P = unscale_state_covariance( S , Pz )

P = S * Pz * S;

end
