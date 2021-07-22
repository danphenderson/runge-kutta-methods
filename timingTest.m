% This script performs a race between our dynamic rk
% methods verse a preallocated rk method
% The test is one three random systems.

clear;

% Define a few functions
a = 0.0001;
b = 1/14;
Pop = 10^4;
F = @(Y)[-a*Y(1)*Y(2); a*Y(1)*Y(2)-b*Y(2)];   

% Redefine functions so we can run them in the matlab solvers
f = @(t,Y)[-a*Y(1)*Y(2); a*Y(1)*Y(2)-b*Y(2)];
options = odeset('RelTol',1e-8,'AbsTol',1e-7);

% Timing 4-5 rk's and ode45..
tic
fprintf('Computation of rk45Dynamic: ')
for i = 1:100
[T1,Y1]=rk45Dynamic(F, 0, 600, [Pop-1; 1], 1e-7);
end
toc

tic
fprintf('Computation of rk45: ')
for i = 1:100
rk45(F, 0, 600, [Pop-1; 1], 1e-7);
end
toc

tic
fprintf('Computation of ode45 ')
for i = 1:100
[t1,y1]=ode45(f, [0 600], [Pop-1; 1],options);
end
toc

fprintf('\n')

% Timing 2-3 rk's and ode23..
tic
fprintf('Computation of rk23Dynamic ')
for i = 1:100
[T2,Y2]=rk23Dynamic(F, 0, 600, [Pop-1; 1], 1e-6);
end
toc

tic
fprintf('Computation of rk23 ')
for i = 1:100
rk23(F, 0, 600, [Pop-1; 1], 1e-6);
end
toc

tic
fprintf('Computation of ode23 ')
for i = 1:100
[t2,y2]=ode23(f, [0 60], [Pop-1; 1],options);
end
toc

% Comparison of step size in the 23 methods
[t1,y1]=ode23(f, [0 60], [Pop-1; 1],options);
[y2,t2]=rk23(F, 0, 60, [Pop-1; 1], 1e-6);
[y3,t3]=rk45(F, 0, 60, [Pop-1; 1], 1e-7);
[t4,y4]=ode45(f, [0 60], [Pop-1; 1],options);

stepChange1 = 0;
stepChange2 = 0;
stepChange3 = 0;
stepChange4 = 0;

% ode23
for i = 2:size(t1, 1)
    stepChange1(i) = t1(i)-t1(i-1);
end
% rk23
for i = 2:size(t2, 2)
    stepChange2(i) = t2(1,i)-t2(1,i-1);
end
% rk45
for i = 2:size(t3, 2)
    stepChange3(i) = t3(1,i)-t3(1,i-1);
end
% ode45
for i = 2:size(t4, 1)
    stepChange4(i) = t4(i)-t4(i-1);
end
figure(1)
semilogy(t1(:,1)',stepChange1, t2(1,:), stepChange2, t3(1,:), stepChange3, t4(:,1)',stepChange4)
legend('ode23','rk23','rk45', 'ode45')
ylabel('Log of Step Change')
xlabel('Time')
title('Plot of Step Change verse Time For Each Method')





