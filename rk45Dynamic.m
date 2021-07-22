function [U,t] = rk45Dynamic(f,t0,tf,u0,tol)
% Solves the initial value problem defined by the input function f,
% where f is an atonomous system of n first-order ODE's having elements 
% of the form u'=f(u), with the initial condition in u0. Both f and u0
% are column vectors of length n. The solution is assumed to exist over 
% the time interval t0 to tf. The relative local truncation error is
% estimated at each step and it must be less than specified tol to move
% forward in time.

% initialize
capacity = (tf-t0)*100;
dt = 0.1;                          % initial step size
i = 1;                              % current solution step                  
U = zeros(size(u0,1), capacity);    % solution matrix
t = zeros(1, capacity);             % points of approximation
U(:,i) = u0;                        % store initial condition
t(i) = t0;                          % store initial starting time

tol = tol/(tf-t0);
while (t(i) < tf) && (i < 1e7)
    % perform numerical approximation
    k1 = f(U(:,i))*dt;
    k2 = f(U(:,i)+0.2*k1)*dt;
    k3 = f(U(:,i)+1/40*(3*k1+9*k2))*dt;
    k4 = f(U(:,i)+44/45*k1-56/15*k2+32/9*k3)*dt;
    k5 = f(U(:,i)+19372/6561*k1-25360/2187*k2+64448/6561*k3-212/729*k4)*dt;
    k6 = f(U(:,i)+9017/3168*k1-355/33*k2+46732/5247*k3+49/176*k4-5103/18656*k5)*dt;
    rk4 = U(:,i)+35/384*k1+500/1113*k3+125/192*k4-2187/6784*k5+11/84*k6;
    k7 = f(rk4)*dt;
    rk5 = U(:,i)+5179/57600*k1+7571/16695*k3+393/640*k4-92097/339200*k5+187/2100*k6+1/40*k7;
    % handle step control
    [U,t,dt,i] = StepControl(rk5,rk4,U,t,i,dt,4,tol);
    % resize solution if necessary
    if capacity == i
        [U,t] = DoubleSize(U,t);
        capacity = capacity*2;
    end
end
U = U(:,1:i);
t = t(:,1:i);