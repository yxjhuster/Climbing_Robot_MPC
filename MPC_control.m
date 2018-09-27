clc;clear;
x_initial = [0;0;0;0];
x_ref = [0.6*pi;0;-0.6*pi;0];
Ts = 2;
N = 40;
uopt = zeros(N,1);
u = [];
options = optimoptions('fmincon','Algorithm','sqp','Display','iter');
max = x_ref(1);
xi = x_initial;
x_h = x_initial;
i=40;
while i>0
N = i;
COSTFUN = @(u) Climbing_Obj(xi,u,Ts,N,x_ref);
CONSFUN = @(u) Climbing_Con(xi,u,Ts,N,max);
uopt = fmincon(COSTFUN,uopt,[],[],[],[],[],[],CONSFUN,options); %need to add more constraints
u = [u uopt(1)];
xk1 = Climbing_DT(xi,uopt(1),Ts);
x_h = [x_h xk1];
xi = xk1(:,end);
i=i-1;
Ts = Ts - 1/40;
end