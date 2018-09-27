function [c, ceq] = Climbing_Con(x,u,Ts,N,max)
c = zeros(N,1);
xk = x;
uk = u(1);
for ct=1:N-1
    xk1 = Climbing_DT(xk,uk,Ts);
    c(ct+1) = xk1(1)-max;
    if ct==1
        c(ct)=x(1);
    end
end
%% No equality constraints
ceq = [];