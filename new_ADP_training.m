clc;clear;
%% Initialization of the system
Ts = 0.1; % sample time
Time = 2; % total duration
N = Time / Ts; % step numbers
xk = [(rand(1)-1/2) *2 * pi; (rand(1)-1/2) * 2 * pi; (rand(1)-1/2) * 2 * pi; (rand(1)-1/2) * 2 * pi];
u = 20 * (rand(1) - 1/2);
phi = @(X) [1; X(1); X(2); X(3); X(4);
                X(1)^2; X(2)^2; X(3)^2; X(4)^2;
                X(1)*X(2); X(1)*X(3); X(1)*X(4); X(2)*X(3); X(2)*X(4); X(3)*X(4); 
                X(1)^3; X(2)^3; X(3)^3; X(4)^2; 
                X(1)^2*X(2); X(1)*X(2)^2; X(1)^2*X(3); X(1)*X(3)^2; X(1)^2*X(4); X(1)*X(4)^2; 
                X(2)^2*X(3); X(2)*X(3)^2; X(2)^2*X(4); X(2)*X(4)^2; 
                X(3)^2*X(4); X(3)*X(4)^2;
                X(1)*X(2)*X(3); X(1)*X(2)*X(4); X(1)*X(2)*X(4); X(2)*X(3)*X(4)
               ];
% phi_cost = @(X, u) [1; X(1); X(2); X(3); X(4); u;
%                 X(1)^2; X(2)^2; X(3)^2; X(4)^2; u^2;
%                 X(1)*X(2); X(1)*X(3); X(1)*X(4); X(1)*u; X(2)*X(3); X(2)*X(4); X(2)*u; X(3)*u; X(3)*X(4); X(4)*u;
%                 X(1)^3; X(2)^3; X(3)^3; X(4)^3; u^3;
%                 X(1)^2*X(2); X(1)*X(2)^2; X(1)^2*X(3); X(1)*X(3)^2; X(1)^2*X(4); X(1)*X(4)^2; X(1)^2*u; X(1)*u^2; 
%                 X(2)^2*X(3); X(2)*X(3)^2; X(2)^2*X(4); X(2)*X(4)^2; X(2)^2*u; X(2)*u^2; 
%                 X(3)^2*X(4); X(3)*X(4)^2; X(3)^2*u; X(3)*u^2; 
%                 X(4)^2*u; X(4)*u^2; 
%                 X(1)*X(2)*X(3); X(1)*X(2)*X(4); X(1)*X(2)*u;  X(1)*X(3)*X(4); X(1)*X(3)*u;X(1)*X(4)*u; X(2)*X(3)*X(4); X(2)*X(3)*u; X(2)*X(4)*u; X(3)*X(4)*u; 
%                ];   
           
           
%% Training Parameters initialization
FinalW_cost = zeros(length(phi(xk)),N); % the parameters for phi, 4 is the size of the state
FinalW_torque = zeros(length(phi(xk)),N);
xr = [0.6*pi;0;-0.6*pi;0];
Q = diag([5000,1,1,1]); % Weight for final state error
R = 20; %Weight scaler for power
NoOfEquations = 100;

RHS_J = zeros(NoOfEquations,1);
LHS_J = zeros(NoOfEquations,length(phi(xk)));

RHS_U = zeros(NoOfEquations,1);
LHS_U = zeros(NoOfEquations,length(phi(xk)));
%% Training
for t = 0 : N-1
    k = N - t;
    states = []
    for i = 1 : NoOfEquations
        xk = [(rand(1)-1/2) *2 * pi; (rand(1)-1/2) * 2 * pi; (rand(1)-1/2) * 2 * pi; (rand(1)-1/2) * 2 * pi];
        states = [states;xk'];
        if k == N
            J_k_t = (xk - xr)'*Q*(xk - xr);
            u = 20 * (rand(1) - 1/2);
            RHS_J(i,:) = J_k_t;
            LHS_J(i,:) = phi(xk)';
%             LHS_J(i,:) = phi_cost(xk, 0)'; %not working because it is not invertible. one possible solution is to use phi_torque instead of phi_cost.
        else
            weight_cost = FinalW_cost(:, k+1);
            weight_torque = FinalW_torque(:, k+1);
            [umax, J_k_t] = torque_optimizer(xk,weight_cost, weight_torque, Ts, R);
            
            RHS_J(i,:) = J_k_t;
            LHS_J(i,:) = phi(xk)';
            
            RHS_U(i,:) = umax;
            LHS_U(i,:) = phi(xk)';
        end
    end
    
    if k == N
%         FinalW_cost(:,k) = (LHS_J'*LHS_J)^-1*LHS_J'*RHS_J;
        FinalW_cost(:,k) = pinv(LHS_J)*RHS_J;
        FinalW_torque(:,k) = zeros(length(phi(xk)),1);
    else
        if det(LHS_J'*LHS_J)==0 || det(LHS_U'*LHS_U)==0
            fprintf('det phi = 0\n');
            break;
        end
%         FinalW_cost(:,k) = (LHS_J'*LHS_J)^-1*LHS_J'*RHS_J;
        FinalW_cost(:,k) = pinv(LHS_J)*RHS_J;
%         FinalW_torque(:,k) = (LHS_U'*LHS_U)^-1*LHS_U'*RHS_U;
        FinalW_cost(:,k) = pinv(LHS_U)*RHS_J;
    end
    FinalW_cost(:,k)'*phi(xk)
    RHS_J(NoOfEquations, :)
    error_cost = norm(FinalW_cost(:,k)'*phi(xk) - RHS_J(NoOfEquations, :))
        
    if isnan(FinalW_cost(:,k))
        fprintf('Training W_cost is diverging...\n');
        diverged = 1;
        break;
    end
    
    if isnan(FinalW_torque(:,k))
        fprintf('Training W_torque is diverging...\n');
        diverged = 1;
        break;
    end
    
end