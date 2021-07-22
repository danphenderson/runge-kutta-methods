function [u,t,dt,i] = StepControl(rkh, rkl, u, t, i, dt, q, tol)
    % estimate absolute and relative local trucation error
    errEstAbs = max(abs(rkh - rkl));
    errEstRel = max(abs(rkh - rkl)./abs(rkh));
    % set error to be the smallest of the two
    err = min(errEstAbs, errEstRel);
    % if true, decrease step size and repeat ith iteration
    if err > tol
        dt = 0.5*dt;
    % else the current step is good and calculate next step
    else
        t(i+1) = t(i)+dt;
        u(:,i+1) = rkh;
        i = i+1;
        dt = 0.9*dt*(tol/err)^(1/q); 
    end
end