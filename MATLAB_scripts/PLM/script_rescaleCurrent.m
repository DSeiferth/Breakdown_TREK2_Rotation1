for mV = [-200 -160:20:-20]
save=1
%fixedEdges = -1.2:0.015:0.2;
fixedEdges = -1.2:0.005:0.2;
%folder = 'rescaled/minScale/';
%folder = 'rescaled/modeScale/fixedEdges/';
folder = 'rescaled/stdScale/fixedEdges/';
%folder = 'rescaled/stdScale/logScale/';
subset = 1e4;
disp(append('Voltage = ',int2str(mV)));
if mV==200
    name = append('../../../SUBSTATES_200/+200.txt');
    subset = 2 * subset;
elseif mV==-200
    name = append('../../../SUBSTATES_200/-200.txt');
    subset = 0.4 * subset;
elseif mV == 0
    name = append('../../../SUBSTATES/zero.txt');
elseif mV < 0
    name = append('../../../SUBSTATES/',int2str(abs(mV)) ,'.txt');
elseif mV > 0
    name = append('../../../SUBSTATES_positive/',int2str(mV) ,'.txt');
    subset = 2 * subset;
end

table = readtable(name);
conv =  table2array(table);
close all

clear table;
% check if channel inverted
if mV<200 && mV >-200
    x_1 = -conv(:,2);
else
    x_1 = conv(:,2);
end
%%%%%%%%%%%%%%%%%%%%% Re-scaling  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%determine local maxima
[N,edges] = histcounts(x_1, 'Normalization','pdf');
[TF1,P] = islocalmax(N, 'MaxNumExtrema', 2);
posMode = edges(TF1)
posOpenMode = posMode(1);
posClosedMode = posMode(2);

minCurrent = min(x_1)
maxCurrent = max(x_1)
rangeCurrent = -minCurrent+maxCurrent
stdCurrent=std(x_1)
stdRange = 3*stdCurrent
%rescaled = (-x_1 + maxCurrent)/rangeCurrent;
%rescaled = (x_1-posClosedMode)/abs(minCurrent);
%rescaled = (x_1-posClosedMode)/abs(posOpenMode-posClosedMode);
rescaled = (x_1-posClosedMode)/abs(stdRange);

%%%%%%%%%%%%   Plot time series %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
step = 5e-6;
t_1 = [0:step:step*(subset-1)];
offset = 1e5;
plot(t_1, rescaled(1+offset:subset+offset))
xlabel('Time [s]')
ylabel('Re-scaled Current')
title(append(int2str(mV), 'mV '))
set(gca,'FontSize',14)
plotname = append(folder, 'timeseries_',int2str(mV) ,'mV');
if save == 1
saveas(gcf,plotname,'epsc');
end

%%%%%%%%%%%%%%%%    plot histograms        %%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
[N,edges] = histcounts(rescaled, fixedEdges,  'Normalization','pdf');
edges = edges(2:end) - (edges(2)-edges(1))/2;
plot(edges, N);
ylabel('pdf')
xlabel('Re-scaled Current')
xlim([-1.2, 0.2])
set(gca, 'YScale', 'log')

%determine local maxima
[TF1,P] = islocalmax(N, 'MaxNumExtrema', 2);
posMode = edges(TF1)

title({[int2str(mV), 'mV']
    ['open mode at ', num2str(posMode(1)) ]
    ['closed mode at ', num2str(posMode(2)) ]
    })
set(gca,'FontSize',14)
plotname = append(folder, 'distribution_',int2str(mV) ,'mV');
if save == 1
saveas(gcf,plotname,'epsc');
fname1 = append(folder, num2str(mV),'mV_N.txt');
fname2 = append(folder, num2str(mV),'mV_edges.txt');
writematrix(N', fname1)
writematrix(edges', fname2)
end

%%%%%%%%%%%%   Estimate noise level   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find region where long in closed state
closedLevel = rescaled > - 0.2;
% https://uk.mathworks.com/matlabcentral/answers/477407-how-to-get-the-longest-consecutive-values-in-a-column-vector-and-the-position-at-which-it-starts
i=reshape(find(diff([0;closedLevel;0])~=0),2,[]);
[lgtmax,jmax]=max(diff(i));
istart=i(1,jmax);
lgtmax % length of the longest sequence of 1s
istart % where it starts

t_c = [0:step:step*(lgtmax-1)];
longestClosedLevel = rescaled(istart:istart+lgtmax-1);
figure(3)
plot(t_c, longestClosedLevel)
title({[int2str(mV), 'mV ']
    ['channel in closed state - noise estimation']})
xlabel('Time [s]')
ylabel('Re-scaled Current')
set(gca,'FontSize',14)
plotname = append(folder, 'noiseLevel_',int2str(mV) ,'mV');
if save == 1
saveas(gcf,plotname,'epsc');
end

figure(4)
[N,edges] = histcounts(longestClosedLevel, 'Normalization','pdf');
edges = edges(2:end) - (edges(2)-edges(1))/2;
plot(edges, N);
ylabel('pdf')
xlabel('Re-scaled Current')
stdClosed = std(longestClosedLevel);
title(append('Noise of closed state, \sigma = ', num2str( stdClosed )))
set(gca,'FontSize',14)
plotname = append(folder, 'noiseLevel_distr_',int2str(mV) ,'mV');
if save == 1
saveas(gcf,plotname,'epsc');
end
%%%%%%%%%%%%   Output  modes and sigma  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save == 1
out_name = append(folder, 'Breakdown_Modes_Sigma.txt');
fileID = fopen(out_name,'a');
fprintf(fileID, '%d, %d, %d, %d\n',mV,posMode(1),posMode(2),stdClosed);
fclose(fileID);
end
end