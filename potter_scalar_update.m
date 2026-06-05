%-------------------------------------------------------------------------%
%    potter_scalar_update
%
%  DESCRIPTION:
%     Potter square-root measurement update for one scalar observation.
%
%  MAE 182 -- Spacecraft Guidance and Navigation
%-------------------------------------------------------------------------%
function [ W , x ] = potter_scalar_update( W , x , h , innov , r )

f     = W' * h';
alpha = f' * f + r;
gam   = 1 / ( alpha + sqrt( alpha * r ) );
b     = W * f;
W     = W - gam * b * f';
x     = x + ( b / alpha ) * innov;

end
