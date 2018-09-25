function J = Climbing_Obj(x, u, Ts, N, xref)
%% Objective function for MPC
% system current states (x): 4x1 vector
% motor torque (u): Nx1 vector
% sample time (Ts): scalar
% number of control sequence (N): scalar
% reference state (xref): 4x51 vector
Q = diag([5000,500,1,500]); %Weight matrix for states error
R = 20; %Weight scaler for power
xk = x;
uk = u(1);
J = 0;
% compute cost for power consumption
for ct=1:N-1
    xk1 = Climbing_DT(xk, uk, Ts);
    J = J + R*(uk*xk1(4)*Ts)^2;
    xk = xk1;
    uk = u(ct+1);
end
% compute cost for state error
xk1 = Climbing_DT(xk, uk, Ts);
J = J + (xk1 - xref)'*Q*(xk1 - xref);
end