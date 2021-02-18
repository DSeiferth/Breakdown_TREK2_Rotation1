close all, clear all
for mV = -100%[-200 -160:20:-20]
disp(append('Voltage = ',int2str(mV)));
folder = 'rescaled/stdScale/fixedEdges/'; 
save = 1;

%%%%%%%%%%%%%%%%%%%%%%   Read Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model_shift = -1;
% stdSclae
fname_N = append(folder,num2str(mV),'mV_N.txt')
fname_edges = append(folder,num2str(mV),'mV_edges.txt')
exp_N = readmatrix(fname_N);
exp_edges = readmatrix(fname_edges);

%determine local maxima
[TF1,P] = islocalmax(exp_N, 'MaxNumExtrema', 2);
posMode = exp_edges(TF1)
d1 = posMode(1)+1
if d1>0.3
    d1=0.3;
end

%%%%%%%%%%%%%%%%%%%%%% Parameters, generate trajectory  %%%%%%%%%%%%%%%%
if mV == -160
    a1=3.1; a2=1.01; a3=1.0005; d1=.17; d2=0.8; d3=0.9995;  %3-states, 1-switching 160mV
elseif mV == -100
    a1=3.1; a2=1.01; a3=1.01; d2=0.75; d3=0.99;  %3-states, 1-switching 100mV
elseif mV == -40
    a1=3.45; a2=1.007; a3=1.05; d2=0.6; d3=0.98;  %3-states, 1-switching 40mV
end
start = 0.99;
sigma = 0.02;

No = 500000;
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
%%%%%%%%%%%%%%%%%%%%%%%%%  PLOT  time series  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
subset = 700;
shift=1e3;
plot(shift:1:shift+subset, noisy_series(shift:1:shift+subset)+model_shift)
hold on
plot(shift:1:shift+subset, timeseries(shift:1:shift+subset)+model_shift)
legend( append('White noise added with \sigma=',num2str(sigma)), 'Iterated map')
xlim([shift, shift+subset])
ylim([-0.1 1.35]+model_shift)
ylabel('Current')
xlabel('Time steps')
title({
    [num2str(mV) 'mV']
    }, 'FontWeight', 'normal')
set(gca,'FontSize',15)
if save == 1
plot_name = append(folder,'model/',num2str(mV),'_mV_time_series');
saveas(gcf,plot_name,'epsc');
end

%%%%%%%%%%%%%%%%%%%%%%    Plot distribution   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
startEdges = (-1.2:0.015:0.2)+1;
[N2,edges2] = histcounts(noisy_series, startEdges, 'Normalization','pdf');
edges2 = edges2(2:end) - (edges2(2)-edges2(1))/2;
plot(edges2+model_shift, N2, '-');
hold on 
[N,edges] = histcounts(timeseries, startEdges, 'Normalization','pdf');
edges = edges(2:end) - (edges(2)-edges(1))/2;
plot(edges+model_shift, N, '-');
hold on
plot(exp_edges, exp_N, 'color', 'black')
legend('Iterated map with noise','Iterated map: no noise','Re-scaled exp. data ', 'Location', 'NorthWest')
ylabel('PDF')
xlabel('Current')
xlim([-1.2, 0.2])
ylim([0, 12])
set(0,'DefaultAxesTitleFontWeight','normal');
title({
    [num2str(mV) 'mV, noise \sigma=',num2str(sigma)]
    ['1st closed interval -1 < x < ' num2str(d1+model_shift)]
    ['2nd closed interval ' num2str(d1+model_shift) ' < x < '  num2str(d2+model_shift)]
    ['open interval ' num2str(d3+model_shift) ' < x < 0']
    ['Slopes: a_1=' num2str(a1) ', a_2=' num2str(a2) ', a_3=' num2str(a3)]
    }, 'FontWeight', 'normal')
set(gca,'FontSize',14)
if save == 1
plot_name = append(folder,'model/',num2str(mV),'mV_amplitude_distr');
saveas(gcf,plot_name,'epsc')
end
%%%%%%%%%%%%%%%%%%%% Plot distribution with log-scale   %%%%%%%%%%%%%%%%%%
figure(3)
[N2,edges2] = histcounts(noisy_series, startEdges, 'Normalization','pdf');
edges2 = edges2(2:end) - (edges2(2)-edges2(1))/2;
plot(edges2+model_shift, N2, '-');
hold on 
[N,edges] = histcounts(timeseries, startEdges, 'Normalization','pdf');
edges = edges(2:end) - (edges(2)-edges(1))/2;
plot(edges+model_shift, N, '-');
hold on
plot(exp_edges, exp_N, 'color', 'black')
legend('Iterated map with noise','Iterated map: no noise','Re-scaled exp. data ', 'Location', 'South')
ylabel('PDF')
xlabel('Current')
xlim([-1.2, 0.2])
set(0,'DefaultAxesTitleFontWeight','normal');
title({
    [num2str(mV) 'mV, noise \sigma=',num2str(sigma)]
    ['1st closed interval -1 < x < ' num2str(d1+model_shift)]
    ['2nd closed interval ' num2str(d1+model_shift) ' < x < '  num2str(d2+model_shift)]
    ['open interval ' num2str(d3+model_shift) ' < x < 0']
    ['Slopes: a_1=' num2str(a1) ', a_2=' num2str(a2) ', a_3=' num2str(a3)]
    }, 'FontWeight', 'normal')
set(gca,'FontSize',14)
set(gca, 'YScale', 'log')
if save == 1
plot_name = append(folder,'model/',num2str(mV),'mV_amplitude_distr_log');
saveas(gcf,plot_name,'epsc')
end
%%%%%%%%%%%%%%%%%%%   Plot map   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)
x = 0:0.001:1;
f_x = zeros(length(x),1);
i=1
for xi=x
    f_x(i) = piecewise3states(xi, a1,a2,a3,d1,d2,d3);
    i=i+1;
end
plot(x+model_shift,f_x+model_shift,'.')
hold on
plot([0,1]+model_shift,[0,1]+model_shift,':', 'color', 'black')
hold on 
plot([d1,d1]+model_shift, [0,1]+model_shift,'--','color', 'black')
hold on 
plot([d2,d2]+model_shift, [0,1]+model_shift,'--','color', 'black')
hold on 
plot([d3-0.005,d3-0.005]+model_shift, [0,1]+model_shift,'--','color', 'black')
ylabel('x_{n+1} = f(x_n)')
xlabel('x_n')
xlim([0, 1]+model_shift)
ylim([0, 1]+model_shift)
set(gca,'FontSize',14)
set(0,'DefaultAxesTitleFontWeight','normal');
title(append('Iterated map for ', num2str(mV), 'mV'), 'FontWeight', 'normal')
if save == 1
plot_name = append(folder,'model/',num2str(mV),'mV_map');
saveas(gcf,plot_name,'epsc')
end
end

