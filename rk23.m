function [U,t] = rk23(f,t0,tf,u0,tol)
% Solves the initial value problem defined by the input function f,
% where f is an atonomous system of n first-order ODE's having elements 
% of the form u'=f(u), with the initial condition in u0. Both f and u0
% are column vectors of length n. The solution is assumed to exist over 
% the time interval t0 to tf. The relative local truncation error is
% estimated at each step and it must be less than specified tol to move
% forward in time.
%
% fixed allocation - approximation stops after 1e4-1 steps

% initialize
dt = 0.01;                           % initial step size
i = 1;                               % current solution step                 
U = zeros(size(u0,1), 1e5);          % solution matrix
t = zeros(1, 1e5);                   % points of approximation
U(:,i) = u0;                         % store initial condition
t(i) = t0;                           % store initial starting time

tol = tol/(tf-t0);
% perform numerical approximation
while (t(i) < tf) && (i < 1e5)
    m1 = f(U(:,i))*dt;
    m2 = f(U(:,i)+0.5*m1)*dt;
    m3 = f(U(:,i)+0.75*m2)*dt;
    rk2 = U(:,i)+2/9*m1+1/3*m2+4/9*m3;
    m4 = f(rk2)*dt;
    rk3 = U(:,i)+7/24*m1+1/4*m2+1/3*m3+1/8*m4;
    % handle step control
    [U,t,dt,i] = StepControl(rk3,rk2,U,t,i,dt,2,tol);
end
U = U(:,1:i);
t = t(:,1:i);