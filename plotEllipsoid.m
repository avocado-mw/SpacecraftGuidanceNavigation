%-------------------------------------------------------------------------%
%    plotEllipsoid
%
%  DESCRIPTION:
%     Plots a 3D covariance ellipsoid given an orthonormal eigenvector
%     matrix R and a vector of semi-axis lengths.
%
%  MAE 182 -- Spacecraft Guidance and Navigation
%-------------------------------------------------------------------------%
function plotEllipsoid( R , semi )

semi = real( semi(:) );
semi = max( semi , 0 );

[x,y,z] = sphere( 24 );

x = x * semi(1);
y = y * semi(2);
z = z * semi(3);

[m,n] = size( x );
C     = ( R * [ x(:) y(:) z(:) ].' ).';

x = reshape( C(:,1) , m , n );
y = reshape( C(:,2) , m , n );
z = reshape( C(:,3) , m , n );

surf( x , y , z , 'FaceAlpha' , 0.35 , 'EdgeColor' , [0.3 0.3 0.3] )
axis equal
grid on

end
