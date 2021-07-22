% Problem: 
% x'(t) = -x + 30e^(-t)cos(30t) + cos(t) + sin(t), x(0) = 0
% which has an exact solution of x(t) = e^(-t)sin(30t) + sin(t)
%
% we encode the problem into an array X by letting
% X(1) = t
% X(2) = x
% then differenting the system gives
% X(1)' = 1
% X(2)' = x' = -X(2) + 30e^(-X(1))cos(30X(1)) + cos(X(1)) + sin(X(1))

clear

% Define ODE system and exact solution
f = @(X)[1; -X(2) + 30*exp(-X(1))*cos(30*X(1)) + cos(X(1)) + sin(X(1))];
F = @(t,X)[1; -X(2) + 30*exp(-X(1))*cos(30*X(1)) + cos(X(1)) + sin(X(1))];
exactSol = @(t)exp(-t).*sin(30.*t) + sin(t);
options = odeset('RelTol',1e-8);


% calculate some inf-norms
[sol1, t1] = rk45(f,0,15,[0; 0],1e-8);
e1 = max(abs(sol1(2,:) - exactSol(t1)));
r1 = norm((sol1(2,:) - exactSol(t1)),2);
fprintf('Euclidean Norm of Error vecotor for rk45: %e\n', r1)
fprintf('Infinity Norm of Error vector for rk45: %e\n\n',e1)


[sol2, t2] = rk45Dynamic(f,0,15,[0; 0],1e-8);
e2 = max(abs(sol2(2,:) - exactSol(t2)));
r2 = norm((sol2(2,:) - exactSol(t2)),2);
fprintf('Euclidean Norm of Error vecotor for rk45Dynamic: %e\n', r2)
fprintf('Infinity Norm of Error vector for rk45Dynamic: %e\n\n',e2)

[T2, Sol2] = ode45(F, [0 15], [0; 0], options);
e6 = max(abs(Sol2(:,2) - exactSol(T2)));
r6 = norm((Sol2(:,2) - exactSol(T2)),2);
fprintf('Euclidean Norm of Error vector for ode45: %e\n', r6)
fprintf('Infinity Norm of Error vector for ode45: %e\n\n', e6)

[sol3, t3] = rk23(f,0,15,[0; 0],1e-8);
e3 = max(abs(sol3(2,:) - exactSol(t3)));
r3 = norm((sol3(2,:) - exactSol(t3)),2);
fprintf('Euclidean Norm of Error vecotor for rk23: %e\n', r3)
fprintf('Infinity Norm of Error vector for rk23:%e\n\n', e3)

[sol4, t4] = rk23Dynamic(f,0,15,[0; 0],1e-8);
e4 = max(abs(sol4(2,:) - exactSol(t4)));
r4 = norm((sol4(2,:) - exactSol(t4)),2);
fprintf('Euclidean Norm of Error vecotor for rk23Dynamic: %e\n', r4)
fprintf('Infinity Norm of Error vector for rk23Dynamic: %e\n\n', e4)


[T1, Sol1] = ode23(F, [0 15], [0; 0], options);
e5 = max(abs(Sol1(:,2) - exactSol(T1)));
r5 = norm((Sol1(:,2) - exactSol(T1)),2);
fprintf('Euclidean Norm of Error vector for ode23: %e\n', r5)
fprintf('Infinity Norm of Error vector for ode23: %e\n\n', e5)


% ode23
for i = 2:size(T1, 1)
    stepChange3(i) = T1(i)-T1(i-1);
end
% ode45
for i = 2:size(T2, 1)
    stepChange4(i) = T2(i)-T2(i-1);
end

% Plot step verse step
for i = 2:size(t2, 2)
    stepChange1(i) = t2(1,i)-t2(1,i-1);
end
figure(1)
semilogy(t2, stepChange1, T2, stepChange4)
title('Step Size of Transient Function Approximation For ode45, rk45.m and rk45Dynamic.m')
ylabel('Log of Step Size')
xlabel('Time')
legend('rk45 and rk45Dynamic', 'ode45')

for i = 2:size(t4, 2)
    stepChange2(i) = t4(1,i)-t4(1,i-1);
end
figure(2)
semilogy(t4,stepChange2, T1, stepChange3)
title('Plot of Step Size of Transient Function Approximation For ode23, rk23.m and rk23Dynamic.m')
ylabel('Log of Step Size')
xlabel('Time')
legend('rk23 and rk23Dynamic', 'ode23')
