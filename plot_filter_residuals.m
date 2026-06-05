%-------------------------------------------------------------------------%
%    plot_filter_residuals
%
%  DESCRIPTION:
%     Plot pre-fit and post-fit residuals versus observation epoch [s].
%     Residuals are grouped by measurement type so the full arc variation
%     is visible even when only one tracking station is used.
%
%  MAE 182 -- Spacecraft Guidance and Navigation
%-------------------------------------------------------------------------%
function plot_filter_residuals( t_obs , pre_fit , post_fit , figPrefix )

figure( 'Name' , [ figPrefix , ' Range Residuals' ] );
plot( t_obs , pre_fit(1,:) , 'r.' , 'DisplayName' , 'pre-fit' ); hold on
plot( t_obs , post_fit(1,:) , 'k.' , 'DisplayName' , 'post-fit' );
grid on; legend( 'Location' , 'best' );
title( [ figPrefix , ' Range Residuals' ] );
xlabel( 'epoch [s]' ); ylabel( 'O - C [m]' );

figure( 'Name' , [ figPrefix , ' Range-Rate Residuals' ] );
plot( t_obs , pre_fit(2,:) , 'r.' , 'DisplayName' , 'pre-fit' ); hold on
plot( t_obs , post_fit(2,:) , 'k.' , 'DisplayName' , 'post-fit' );
grid on; legend( 'Location' , 'best' );
title( [ figPrefix , ' Range-Rate Residuals' ] );
xlabel( 'epoch [s]' ); ylabel( 'O - C [m/s]' );

end
