clc;clear;
%% Initialization of the system
Ts = 0.01; % sample time
Time = 5; % total duration
N = Time / Ts; % step numbers
xk = [(rand(1)-1/2) *2 * pi; (rand(1)-1/2) * 2 * pi; (rand(1)-1/2) * 2 * pi; (rand(1)-1/2) * 2 * pi];
u = 20 * (rand(1) - 1/2);
% Up to the Third order % !!!!!!!! I think I need to add u inside the phi.
% phi = @(X) [1; X(1); X(2); X(3); X(4);
%                 X(1)^2; X(2)^2; X(3)^2; X(4)^2;
%                 X(1)*X(2); X(1)*X(3); X(1)*X(4); X(2)*X(3); X(2)*X(4); X(3)*X(4); 
%                 X(1)^3; X(2)^3; X(3)^3; X(4)^2; 
%                 X(1)^2*X(2); X(1)*X(2)^2; X(1)^2*X(3); X(1)*X(3)^2; X(1)^2*X(4); X(1)*X(4)^2; 
%                 X(2)^2*X(3); X(2)*X(3)^2; X(2)^2*X(4); X(2)*X(4)^2; 
%                 X(3)^2*X(4); X(3)*X(4)^2;
%                 X(1)*X(2)*X(3); X(1)*X(2)*X(4); X(1)*X(2)*X(4); X(2)*X(3)*X(4)
%                ];
phi = @(X, u) [1; X(1); X(2); X(3); X(4); u;
                X(1)^2; X(2)^2; X(3)^2; X(4)^2; u^2;
                X(1)*X(2); X(1)*X(3); X(1)*X(4); X(1)*u; X(2)*X(3); X(2)*X(4); X(2)*u; X(3)*u; X(3)*X(4); X(4)*u;
                X(1)^3; X(2)^3; X(3)^3; X(4)^3; u^3;
                X(1)^2*X(2); X(1)*X(2)^2; X(1)^2*X(3); X(1)*X(3)^2; X(1)^2*X(4); X(1)*X(4)^2; X(1)^2*u; X(1)*u^2; 
                X(2)^2*X(3); X(2)*X(3)^2; X(2)^2*X(4); X(2)*X(4)^2; X(2)^2*u; X(2)*u^2; 
                X(3)^2*X(4); X(3)*X(4)^2; X(3)^2*u; X(3)*u^2; 
                X(4)^2*u; X(4)*u^2; 
                X(1)*X(2)*X(3); X(1)*X(2)*X(4); X(1)*X(2)*u;  X(1)*X(3)*X(4); X(1)*X(3)*u;X(1)*X(4)*u; X(2)*X(3)*X(4); X(2)*X(3)*u; X(2)*X(4)*u; X(3)*X(4)*u; 
               ];           
%% Training
FinalW = zeros(length(phi(xk,u)),N); % the parameters for phi, 4 is the size of the state
xr = [0.6*pi;0;-0.6*pi;0];
Q = diag([5000,500,1,500]); % Weight for final state error
R = 20; %Weight scaler for power
NoOfEquations = 1000;
RHS_J = zeros(NoOfEquations,1);%we can take it as y in Ab = y and least square problem is going to find b
LHS_J = zeros(NoOfEquations,length(phi(xk, u))); %we can take it as A in Ab = y and least square problem is going to find b

for t = 0:N-1
    k = N - t;
    for i = 1:NoOfEquations
        xk = [(rand(1)-1/2) *2 * pi; (rand(1)-1/2) * 2 * pi; (rand(1)-1/2) * 2 * pi; (rand(1)-1/2) * 2 * pi];
        u = 10 * (rand(1) - 1/2);
        if k == N
            J_k_t = (xk - xr)'*Q*(xk - xr);
            RHS_J(i,:) = J_k_t;
            LHS_J(i,:) = phi(xk, u)';
        else
            xk_plus_1 = Climbing_DT(xk,u,Ts);
            u_plus_1 = 10 * (rand(1) - 1/2);
            J_k_plus_1 = FinalW(:,k+1)' * phi(xk_plus_1, u_plus_1); %how to find the u_plus_1
            J_k_t = J_k_plus_1 + R * (u * xk(4) * Ts)^2;
            RHS_J(i, :) = J_k_t;
            LHS_J(i, :) = phi(xk, u)';
        end
    end
    
    if det(LHS_J'*LHS_J)==0
        fprintf('det phi = 0\n');
        break;
    end
    
    FinalW(:,k) = (LHS_J'*LHS_J)^-1*LHS_J'*RHS_J;
        
    if isnan(FinalW(:,k))
        fprintf('Training W is diverging...\n');
        diverged = 1;
        break;
    end
end

%the final cost is given by the first column of FinalW times phi.
plot([1:N], FinalW(:,1:N))
xlabel('Time Steps'); ylabel('Weight Elements'); grid on