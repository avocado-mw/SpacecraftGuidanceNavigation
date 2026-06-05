%-------------------------------------------------------------------------%
%    build_P0_bar
%
%  DESCRIPTION:
%     Builds the 18x18 a priori covariance for the MAE 182 estimator.
%     Station 101 (s1) is fixed.  Stations 337 and 394 (s2, s3) share the
%     same initial variance supplied by stationVar_s2s3.
%
%  MAE 182 -- Spacecraft Guidance and Navigation
%-------------------------------------------------------------------------%
function P0_bar = build_P0_bar( stationVar_s2s3 )

if nargin < 1
    stationVar_s2s3 = 1e6;
end

P0_bar = diag([ ...
    1e6 , 1e6 , 1e6 , ...
    1e6 , 1e6 , 1e6 , ...
    1e20, ...
    1e6 , ...
    1e6 , ...
    1e-10, 1e-10, 1e-10, ...
    stationVar_s2s3 , stationVar_s2s3 , stationVar_s2s3 , ...
    stationVar_s2s3 , stationVar_s2s3 , stationVar_s2s3 ]);

end
