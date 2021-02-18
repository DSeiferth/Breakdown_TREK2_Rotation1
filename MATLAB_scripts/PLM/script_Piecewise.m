clear all, close all, clc
%a1=2.19; a2=1.01; d1=0.2; d2=0.95;  % 2-state ala Liebovitch 1991
a1=3.1; a2=1.01; a3=1.01; d1=0.31; d2=0.7,d3=0.99;  %3-states, 1-switching 
start = 0.99;
sigma = 0.02;

folder = '';
save = 0;  %%%% save plots?
out_name = append(folder, 'paramteres.txt');
fileID = fopen(out_name,'a');
%fprintf(fileID, '%d, %d, %d, %d\n',a1,a2,d1,d2);
fclose(fileID);
No = 100000;
timeseries = zeros(No,1);
for n=1:1:No
    %x2 = piecewise(start, a1,a2,d1,d2) %+ normrnd(0,sigma,1,1);
    x2 = piecewise3states(start, a1,a2,a3,d1,d2,d3);
    timeseries(n) = x2;
    start = x2;
end
%disp(timeseries)
noise = normrnd(0,sigma,No,1);
noisy_series = timeseries + noise;
%%%%%%%%%%%%%%%%%%%%%%%%%  PLOT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
subset = 300;
shift=1e3;
plot(shift:1:shift+subset, noisy_series(shift:1:shift+subset))
hold on
plot(shift:1:shift+subset, timeseries(shift:1:shift+subset))
legend( append('White noise added with \sigma=',num2str(sigma)), 'Iterated map')
ylim([-0.1 1.35])
ylabel('Current')
xlabel('Time steps')
set(gca,'FontSize',15)
if save == 1
plot_name = append(folder,'timeseries');
saveas(gcf,plot_name,'epsc');
end

%%%%%%%%%%%%%%%%%%%%%%   Read Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname_N = 'rescaled/modeScale/fixedEdges/-100mV_N.txt'
fname_edges = 'rescaled/modeScale/fixedEdges/-100mV_edges.txt'
exp_N = readmatrix(fname_N);
exp_edges = readmatrix(fname_edges);
% stdSclae
fname_N = 'rescaled/stdScale/fixedEdges/-100mV_N.txt'
fname_edges = 'rescaled/stdScale/fixedEdges/-100mV_edges.txt'
exp_N2 = readmatrix(fname_N);
exp_edges2 = readmatrix(fname_edges);

model_shift = -1;
%%%%%%%%%%%%%%%%%%%%%% Plot distribution  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
startEdges = (-1.5:0.015:0.3)+1;
[N2,edges2] = histcounts(noisy_series, startEdges, 'Normalization','pdf');
edges2 = edges2(2:end) - (edges2(2)-edges2(1))/2;
plot(edges2+model_shift, N2, '-');
hold on 
[N,edges] = histcounts(timeseries, startEdges, 'Normalization','pdf');
edges = edges(2:end) - (edges(2)-edges(1))/2;
%plot(edges+model_shift, N, '-');
hold on
%plot(exp_edges, exp_N)
hold on
plot(exp_edges2, exp_N2)
%legend('Iterated map','exp. data: modeScale ', ...
%    'exp. data: stdScale ', 'Location', 'NorthWest')
ylabel('PDF')
xlabel('Current')
set(0,'DefaultAxesTitleFontWeight','normal');
title({
    ['Open interval 0 < x < ' num2str(d1)]
    ['Closed interval ' num2str(d2) ' < x < 1']
    ['Slopes: a_1=' num2str(a1) ', a_2=' num2str(a2)]
    })
set(gca,'FontSize',14)
if save == 1
plot_name = append(folder,'amplitude_distr');
saveas(gcf,plot_name,'epsc')
end

%%%%%%%%%%%%%% Cost function   %%%%%%%%%%%%%%%%%%%%%
costfunc = sum((N2'-exp_N).^2);