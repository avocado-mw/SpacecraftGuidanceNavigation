%-------------------------------------------------------------------------%
%    propagate_full_arc
%
%  DESCRIPTION:
%     Integrates the reference trajectory and 9x9 STM over the full project
%     time interval, then resamples the results every 20 seconds.
%
%  INPUTS:
%     S      - structure from project_setup
%     X0_ode - ODE initial condition [state; reshape(Phi,81,1)]
%
%  OUTPUTS:
%     t20    - resampled times [s]
%     x20    - resampled full stored state [18 x N]
%     Phi20  - resampled full 18x18 STM history
%
%  MAE 182 -- Spacecraft Guidance and Navigation
%-------------------------------------------------------------------------%
function [ t20 , x20 , Phi20 ] = propagate_full_arc( S , X0_ode )

options = odeset( 'RelTol' , S.tol , 'AbsTol' , S.tol );
[ ~ , X ] = ode45( @project_dyn , S.timeFull , X0_ode , options );

x_temp   = X(:,1:18).';
Phi_temp = X(:,19:end).';
Phi_temp = reshape( Phi_temp , 9 , 9 , numel( S.timeFull ) );

nStore = round( numel( S.timeFull ) / S.rsmpl ) + 1;
t20    = zeros( nStore , 1 );
x20    = zeros( 18 , nStore );
Phi20  = zeros( 18 , 18 , nStore );

for k = 1:nStore
    idx      = S.rsmpl * ( k - 1 ) + 1;
    t20(k)   = S.timeFull(idx);
    x20(:,k) = x_temp(:,idx);
    Phi20(:,:,k) = [ Phi_temp(:,:,idx) , zeros( 9 , 9 ) ; ...
                     zeros( 9 , 9 )    , eye( 9 , 9 )   ];
end

end
