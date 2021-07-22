function [U,T] = DoubleSize(u,t)
% Used to dynamically grow solution size in dynamic runge-kutta methods
% that implements an adaptive time step.
U = zeros(size(u,1),2*size(u,2));
T = zeros(1,2*size(t,2));
U(1:size(u,1), 1:size(u,2)) = u;
T(1,1:size(t,2)) = t;
end