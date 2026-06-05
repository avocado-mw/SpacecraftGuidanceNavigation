%-------------------------------------------------------------------------%
%    propagate_step
%
%  DESCRIPTION:
%     Integrates the reference trajectory and 9x9 STM from ti to tf and
%     returns the state and full 18x18 STM at tf.
%
%  INPUTS:
%     S      - structure from project_setup
%     X0_ode - ODE initial condition [state; reshape(Phi,81,1)] at ti
%     ti     - initial time [s]
%     tf     - final time [s]
%
%  OUTPUTS:
%     x_tf   - stored state at tf [18x1]
%     Phi_tf - full 18x18 STM from ti to tf
%
%  MAE 182 -- Spacecraft Guidance and Navigation
%-------------------------------------------------------------------------%
function [ x_tf , Phi_tf ] = propagate_step( S , X0_ode , ti , tf )

if abs( tf - ti ) < eps
    x_tf   = X0_ode(1:18);
    Phi_tf = eye( 18 );
    return
end

options = odeset( 'RelTol' , S.tol , 'AbsTol' , S.tol );
time    = ti:S.tstep:tf;
if time(end) ~= tf
    time = [ time , tf ];
end

[ ~ , X ] = ode45( @project_dyn , time , X0_ode , options );

x_tf     = X(end,1:18).';
Phi_temp = reshape( X(end,19:end) , 9 , 9 );
Phi_tf   = [ Phi_temp , zeros( 9 , 9 ) ; ...
             zeros( 9 , 9 ) , eye( 9 , 9 ) ];

end
