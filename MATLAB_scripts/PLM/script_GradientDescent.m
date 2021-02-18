close all
%%%%%%%%%%%%%%%%%%%%%%   Read Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mV = -100;
folder = 'rescaled/modeScale/fixedEdges/';
%folder = 'rescaled/stdScale/fixedEdges/';
fname_N = append(folder, num2str(mV), 'mV_N.txt');
fname_edges = append(folder, num2str(mV), 'mV_edges.txt');
exp_N = readmatrix(fname_N);
exp_edges = readmatrix(fname_edges);
model_shift = -1; %for plotting exp and model distribution

a1=2.19; a2=1.01; d1=0.26; d2=0.99;
startParam = [a1,a2,d1,d2];
params1 = startParam;
alpha = 0.00001;
epsilon = 0.0001;
NoIteration = 20;
costArray = zeros(NoIteration, 1);
for i=1:1:NoIteration
    [cost_J, dJ] = Grad_Step(exp_N,params1, epsilon);
    costArray(i) = cost_J;
    params1 = params1 + alpha*dJ;
    if params1(2) < 1
        params1(2) = 1.01;
    end
    if params1(4) > 1
        params1(4) = 0.995;
    end
end
disp(costArray)
disp(params1)
disp(startParam)

%%%%%%%%%%%%  Plot distributions  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
[N2,edges2] = model_distr(startParam(1),startParam(2),startParam(3),startParam(4));
plot(edges2+model_shift, N2, '-');
hold on 
[N2,edges2] = model_distr(params1(1),params1(2),params1(3),params1(4));
plot(edges2+model_shift, N2, '-');
hold on 
plot(exp_edges, exp_N)
legend('initial values','after grad. des.', 'exp. data','Location', 'NorthWest')
ylabel('PDF')
xlabel('Current')
set(0,'DefaultAxesTitleFontWeight','normal');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function J = cost(exp_N, a1, a2, d1, d2)
[model_N, model_edges] = model_distr(a1,a2,d1,d2);
norm=length(exp_N);
J = 1/norm * sum((model_N-exp_N).^2);
end

function [J, dJ] = Grad_Step(exp_N, params, epsilon)
a1 = params(1);
a2 = params(2);
d1 = params(3);
d2 = params(4);
% cost function
J = cost(exp_N,a1,a2,d1,d2);
%derivative
%a1
J_a1_pe = cost(exp_N,a1+epsilon,a2,d1,d2);
J_a1_me = cost(exp_N,a1-epsilon,a2,d1,d2);
dJ_a1 = (J_a1_pe-J_a1_me)/(2*epsilon) ;
%a2
J_a2_pe = cost(exp_N,a1,a2+epsilon,d1,d2);
J_a2_me = cost(exp_N,a1,a2-epsilon,d1,d2);
dJ_a2 = (J_a2_pe-J_a2_me)/(2*epsilon) ;
%d1
J_d1_pe = cost(exp_N,a1,a2,d1+epsilon,d2);
J_d1_me = cost(exp_N,a1,a2,d1+epsilon,d2);
dJ_d1 = (J_d1_pe-J_d1_me)/(2*epsilon) ;
%d2
J_d2_pe = cost(exp_N,a1,a2,d1,d2+epsilon);
J_d2_me = cost(exp_N,a1,a2,d1,d2+epsilon);
dJ_d2 = (J_d2_pe-J_d2_me)/(2*epsilon) ;

dJ = [dJ_a1, dJ_a2, dJ_d1, dJ_d2];
end

function [model_N, model_edges] = model_distr(a1,a2,d1,d2)
sigma = 0.02;
No = 100000;
start = d1/2;
timeseries = zeros(No,1);
for n=1:1:No
    x2 = piecewise(start, a1,a2,d1,d2);
    timeseries(n) = x2;
    start = x2;
end
noise = normrnd(0,sigma,No,1);
noisy_series = timeseries + noise;
startEdges = (-1.5:0.015:0.3)+1;
[N2,edges2] = histcounts(noisy_series, startEdges, 'Normalization','pdf');
model_edges = edges2(2:end) - (edges2(2)-edges2(1))/2;
model_N = N2';
end