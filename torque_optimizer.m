function [umax, J] = torque_optimizer(xk,weight_cost, weight_torque, Ts, R)
max = 10;
min = -10;
umax = 0; %initialize umax
% options = optimoptions('fmincon','Algorithm','sqp','Display','iter');
options = optimoptions('fmincon','Algorithm','sqp');
COSTFUN = @(u) objective(xk, u, Ts, weight_cost, weight_torque, R);
CONSFUN = @(u) con(u, max, min);
[umax, J] = fmincon(COSTFUN,umax,[],[],[],[],[],[],CONSFUN,options); 
end


% function value = phi_cost(X,u)
% value = [1; X(1); X(2); X(3); X(4); u;
%                 X(1)^2; X(2)^2; X(3)^2; X(4)^2; u^2;
%                 X(1)*X(2); X(1)*X(3); X(1)*X(4); X(1)*u; X(2)*X(3); X(2)*X(4); X(2)*u; X(3)*u; X(3)*X(4); X(4)*u;
%                 X(1)^3; X(2)^3; X(3)^3; X(4)^3; u^3;
%                 X(1)^2*X(2); X(1)*X(2)^2; X(1)^2*X(3); X(1)*X(3)^2; X(1)^2*X(4); X(1)*X(4)^2; X(1)^2*u; X(1)*u^2; 
%                 X(2)^2*X(3); X(2)*X(3)^2; X(2)^2*X(4); X(2)*X(4)^2; X(2)^2*u; X(2)*u^2; 
%                 X(3)^2*X(4); X(3)*X(4)^2; X(3)^2*u; X(3)*u^2; 
%                 X(4)^2*u; X(4)*u^2; 
%                 X(1)*X(2)*X(3); X(1)*X(2)*X(4); X(1)*X(2)*u;  X(1)*X(3)*X(4); X(1)*X(3)*u;X(1)*X(4)*u; X(2)*X(3)*X(4); X(2)*X(3)*u; X(2)*X(4)*u; X(3)*X(4)*u; 
%                ];
% end


function value = phi(X)
value =  [1; X(1); X(2); X(3); X(4);
                X(1)^2; X(2)^2; X(3)^2; X(4)^2;
                X(1)*X(2); X(1)*X(3); X(1)*X(4); X(2)*X(3); X(2)*X(4); X(3)*X(4); 
                X(1)^3; X(2)^3; X(3)^3; X(4)^2; 
                X(1)^2*X(2); X(1)*X(2)^2; X(1)^2*X(3); X(1)*X(3)^2; X(1)^2*X(4); X(1)*X(4)^2; 
                X(2)^2*X(3); X(2)*X(3)^2; X(2)^2*X(4); X(2)*X(4)^2; 
                X(3)^2*X(4); X(3)*X(4)^2;
                X(1)*X(2)*X(3); X(1)*X(2)*X(4); X(1)*X(2)*X(4); X(2)*X(3)*X(4)
               ];
end

function J = objective(x, u, Ts, weight_cost, weight_torque, R)
x_plus_1 = Climbing_DT(x, u, Ts);
umax_plus_1 = weight_torque' * phi(x_plus_1);
J_plus_1 = weight_cost' * phi(x_plus_1);
J = J_plus_1 + R*(u*x(4)*Ts)^2;
end


function [c,ceq] = con(u, max, min)
c = zeros(2, 1);
% min <= u <= max
c(1) = u - max;
c(2) = min - u;
ceq = [];
end
