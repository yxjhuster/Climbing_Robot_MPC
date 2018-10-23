close all, clear all, clc, format compact
% number of samples of each class
K = 100;
% define 4 clusters of input data
q = .6; % offset of classes
A = [rand(1,K)-q; rand(1,K)+q];
B = [rand(1,K)+q; rand(1,K)+q];
C = [rand(1,K)+q; rand(1,K)-q];
D = [rand(1,K)-q; rand(1,K)-q];
% plot clusters
figure(1)
plot(A(1,:),A(2,:),'k+')
hold on
grid on
plot(B(1,:),B(2,:),'b*')
plot(C(1,:),C(2,:),'kx')
plot(D(1,:),D(2,:),'bd')
% text labels for clusters
text(.5-q,.5+2*q,'Class A')
text(.5+q,.5+2*q,'Class B')
text(.5+q,.5-2*q,'Class C')
text(.5-q,.5-2*q,'Class D')

% coding (+1/-1) of 4 separate classes
a = [-1 -1 -1 +1]';
b = [-1 -1 +1 -1]';
d = [-1 +1 -1 -1]';
c = [+1 -1 -1 -1]';

% define inputs (combine samples from all four classes)
P = [A B C D];
% define targets
T = [repmat(a,1,length(A)) repmat(b,1,length(B)) ...
 repmat(c,1,length(C)) repmat(d,1,length(D)) ];


net = feedforwardnet([4 3]);
% train net
net.divideParam.trainRatio = 0.9; % training set [%]
net.divideParam.valRatio = 0.05; % validation set [%]
net.divideParam.testRatio = 0.05; % test set [%]
% train a neural network
[net,tr,Y,E] = train(net,P,T);

view(net)

% evaluate performance: decoding network response
[m,i] = max(T); % target class
[m,j] = max(Y); % predicted class
N = length(Y); % number of all samples
k = 0; % number of missclassified samples
if find(i-j), % if there exist missclassified samples
 k = length(find(i-j)); % get a number of missclassified samples
end
fprintf('Correct classified samples: %.1f%% samples\n', 100*(N-k)/N)
% plot network output
figure;
subplot(211)
plot(T')
title('Targets')
ylim([-2 2])
grid on
subplot(212)
plot(Y')
title('Network response')
xlabel('# sample')
ylim([-2 2])
%% 
inputs = [1:6]'; % input vector (6-dimensional pattern)
outputs = [1 2]'; % corresponding target output vector
net = network( ...
1, ...
3, ...
[1;1;1], ...
[1;0;0], ...
[0 0 0; 1, 0, 0; 0,1 0], ...
[0,0,1] ...
);
% view(net)
net.layers{1}.size = 5;
net.layers{2}.size = 5;
net.layers{3}.size = 5;
% view(net)
net = configure(net,inputs,outputs);
initial_output = net(inputs)
net.trainFcn = 'trainlm';
net.performFcn = 'mse';
net.trainParam.showWindow = false;
[net, tr] = train(net,inputs,outputs);
final_output = net(inputs)

