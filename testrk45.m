%% Test Case 2
% SIR Model
% S(t)' = -aS(t)I(t)
% I(t)' = aS(t)I(t) - bI(t)
% R(t)' = bI(t)
% S(t) + I(t) + R(t) = Pop
% we encode this problem into an array Y by letting
% Y(1) = S(t)
% Y(2) = I(t)
% Y(3) = R(t)
% differenting the above system gives
% Y(1)' = S(t)' = -a*Y(1)*Y(2)
% Y(2)' = I(t)' = a*Y(1)*Y(2) - b*Y(2)
% Y(3)' = R(t)' = b*Y(2)
% Note R(t) = Pop - S(t) - I(t) and Y(1)' and Y(2)' are independent
% of Y(3) so we can ignore computing this solution in our rk23 solver.
% Let our parameters be a = 0.0001, b = 1/14, and Pop = 10^4.
% And intially have one person infected which leads to initial conditions:
%       S(0) = Pop - 1 and I(0) = 1
% Then we will solve for t in [0,60]
a = 0.0001;
b = 1/14;
Pop = 10^4;
f2 = @(Y)[-a*Y(1)*Y(2); a*Y(1)*Y(2)-b*Y(2)];
[sol2, t2] = rk45(f2, 0, 60, [Pop-1; 4], 1e-3);
figure(2)
plot(t2, sol2(1,:), t2, sol2(2,:), t2, Pop - sol2(1,:) - sol2(2,:));
