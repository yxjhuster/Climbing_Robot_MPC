function xk1 = Climbing_DT(xk, uk, Ts)
%% Discrete-time dynamic model of climbing robot
% 4 states at k (xk): 
%   mainbody angle (t1)
%   mainbody angular velocity (t1_dot): when positive, moves anti-clockwisely
%   tail angle (t2)
%   tail angular velocity (t2_dot): when positive, moves anti-clockwisely
%
% 1 inputs: (uk)
%   torque (tau): when positive, torque push the tail moves anti-clockwisely 
%
% Sample time (Ts):
%   Sample time for the discretization.
%   In this case, to achieve high accuracy, a further sample rate will
%   added in discretization
%
% 4 states at k+1 (xk1):
%   Same as defined before
%% Parameters
N = 10; % sample rate to achieve high accuracy
delta = Ts / N; % sample time for updating
xk1 = xk; % Initialize the x(k+1)
%% Loop for updating x(k+1)
for ct = 1:N
    xk1 = xk1 + delta * Climbing_CT(xk1, uk);
end
% xk1 = xk1 + Ts * Climbing_DT(xk1, uk); %If you think calculation is slow,
% remove N, and try this command instead of using loop.
end