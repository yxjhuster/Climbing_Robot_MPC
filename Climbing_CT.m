function dxdt = Climbing_CT(x,u)
%% Continuous-time nonlinear dynamic model of climbing robot
% 4 states (x): 
%   mainbody angle (t1)
%   mainbody angular velocity (t1_dot): when positive, moves anti-clockwisely
%   tail angle (t2)
%   tail angular velocity (t2_dot): when positive, moves anti-clockwisely
% 
% 1 inputs: (u)
%   torque (tau): when positive, torque push the tail moves anti-clockwisely 
%
% 4 outputs: (y)
%   same as states (i.e. all the states are measureable)
%   
% dxdt is the derivative of the states.
dxdt = zeros(4,1);
%% Parameters
m1 = 1;  % body mass
m2 = 1;  % tail mass
g = 9.81;   % gravity of earth
l1 = 0.5; % joint length & length from fixed joint to gravity centre of bodym2 * l2^2
l2 = 0.5; % length from joint to the gravity centre of tail
%% Extract parameters from x and u
% x
t1 = x(1);
t1_dot = x(2);
t2 = x(3);
t2_dot = x(4);
% u
tau = u;
%% build the dynamic model with mass matrix M, coriolis matrix C and gravitational matrix N
M = zeros(2,2);
C = zeros(2,2);
N = zeros(2,1);
% assign values for M matrix
M(1,1) = (m1 + m2)*l1^2 + m2 * l2^2 + 2 * m2 * l1 * l2 * cos(t2);
M(1,2) = m2 * l2^2 + m2*l1*l2*cos(t2);
M(2,1) = m2 * l2^2 + m2*l1*l2*cos(t2);
M(2,2) = m2 * l2^2;
% assign values for C matrix
C(1,1) = -2*m2*l1*l2*sin(t2)*t2_dot;
C(1,2) = -m2*l1*l2*sin(t2) * t2_dot;
C(2,1) = m2*l1*l2*sin(t1)*t1_dot;
C(1,1) = 0;
% assign values for N matrix
N(1) = (m1 + m2) * l1 * g * sin(t1) + m2*g*l2*sin(t1+t2);
N(2) = m2*g*l2*sin(t1+t2);
%% Compute the derivatives of the states
% We know that:
%             input = M * state_dot_dot + C * state_dot + N
%It is easy to get:
%             state_dot_dot = inv(M) * (input - N - C * state_dot)
input = [0 ; tau];
state_dot = [t1_dot; t2_dot];
state_dot_dot = inv(M) * (input - N - C* state_dot);
%% assembly all the results
dxdt(1) = t1_dot;
%dxdt(2) = state_dot_dot(1);
dxdt(2) = (t1_dot^2*l1^2*l2*m2*sin(2.0*t2) - 2.0*l1*tau*cos(t2) - 2.0*l2*tau + g*l1*l2*m2*sin(t1 + 2.0*t2) + 2.0*t1_dot^2*l1*l2^2*m2*sin(t2) + 2.0*t2_dot^2*l1*l2^2*m2*sin(t2) - 2.0*g*l1*l2*m1*sin(t1) - 1.0*g*l1*l2*m2*sin(t1) + 4.0*t1_dot*t2_dot*l1*l2^2*m2*sin(t2))/(l1^2*l2*(2.0*m1 + m2 - 1.0*m2*cos(2.0*t2)));
dxdt(3) = t2_dot;
%dxdt(4) = state_dot_dot(2);
dxdt(4) = (tau*(l1^2*m1 + l1^2*m2 + l2^2*m2 + 2.0*l1*l2*m2*cos(t2)) - l1*l2*m2*(t1_dot^2*l1^2*m1*sin(t2) + t1_dot^2*l1^2*m2*sin(t2) + t1_dot^2*l2^2*m2*sin(t2) + t2_dot^2*l2^2*m2*sin(t2) - 1.0*g*l2*m1*sin(t1) - 1.0*g*l2*m2*sin(t1) + 2.0*t1_dot*t2_dot*l2^2*m2*sin(t2) + g*l1*m1*cos(t1)*sin(t2) + g*l1*m2*cos(t1)*sin(t2) + g*l2*m2*cos(t2)^2*sin(t1) + t1_dot^2*l1*l2*m2*sin(2.0*t2) + 0.5*t2_dot^2*l1*l2*m2*sin(2.0*t2) + g*l2*m2*cos(t1)*cos(t2)*sin(t2) + t1_dot*t2_dot*l1*l2*m2*sin(2.0*t2)))/(l1^2*l2^2*m2*(m1 + m2 - 1.0*m2*cos(t2)^2));
end
