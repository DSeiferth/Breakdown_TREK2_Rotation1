clear all, close all, clc
for mV = -200%[-200 -160:20:-20]
save=1
fixedEdges = -70:0.25:10;
folder = 'fixedEdges/';

subset = 1e4;
disp(append('Voltage = ',int2str(mV)));
if mV==200
    name = append('../SUBSTATES_200/+200.txt');
    subset = 2 * subset;
elseif mV==-200
    name = append('../SUBSTATES_200/-200.txt');
    subset = 0.4 * subset;
elseif mV == 0
    name = append('../SUBSTATES/zero.txt');
elseif mV < 0
    name = append('../SUBSTATES/',int2str(abs(mV)) ,'.txt');
elseif mV > 0
    name = append('../SUBSTATES_positive/',int2str(mV) ,'.txt');
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

%%%%%%%%%%%%   Plot time series %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
step = 5e-6;
t_1 = [0:step:step*(subset-1)];
offset = 1e5;
plot(t_1, x_1(1+offset:subset+offset))
xlabel('Time [s]')
ylabel('Current [pA]')
title(append(int2str(mV), 'mV '))
set(gca,'FontSize',14)
plotname = append(folder, 'timeseries_',int2str(mV) ,'mV');
if save == 1
saveas(gcf,plotname,'epsc');
end

%%%%%%%%%%%%%%%%    plot histograms        %%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
[N,edges] = histcounts(x_1, fixedEdges,  'Normalization','pdf');
edges = edges(2:end) - (edges(2)-edges(1))/2;
plot(edges, N);
ylabel('pdf')
ylabel('Current [pA]')
%set(gca, 'YScale', 'log')

% open probability
A = trapz(edges,N)
n = floor(length(edges)/2);
if mV < 0
    close = trapz(edges(n:end),N(n:end))
    open = trapz(edges(1:n),N(1:n))
else
    open = trapz(edges(n:end),N(n:end))
    close = trapz(edges(1:n),N(1:n))
end

%determine local maxima
[TF1,P] = islocalmax(N, 'MaxNumExtrema', 2);
posMode = edges(TF1)
if mV < 0
    posOpenMode = posMode(1);
    posClosedMode = posMode(2);
else
    posOpenMode = posMode(2);
    posClosedMode = posMode(1);
end
title({[int2str(mV), 'mV; P_o = ', num2str(open)]
    ['open mode at ', num2str(posOpenMode) ]
    ['closed mode at ', num2str(posClosedMode) ]
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
disp('posClosedMode * 5')
disp(posClosedMode*5)
if mV < 0
    closedLevel = x_1 > posClosedMode * 5;
else
   closedLevel = x_1 <  posClosedMode * 5;
end
% https://uk.mathworks.com/matlabcentral/answers/477407-how-to-get-the-longest-consecutive-values-in-a-column-vector-and-the-position-at-which-it-starts
i=reshape(find(diff([0;closedLevel;0])~=0),2,[]);
[lgtmax,jmax]=max(diff(i));
istart=i(1,jmax);
lgtmax % length of the longest sequence of 1s
istart % where it starts

t_c = [0:step:step*(lgtmax-1)];
longestClosedLevel = x_1(istart:istart+lgtmax-1);
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
out_name = append(folder, 'mV_Peak1_Peak2_stdClosed_Po.txt');
fileID = fopen(out_name,'a');
fprintf(fileID, '%d, %d, %d, %d, %d\n',mV,posMode(1),posMode(2),stdClosed, open);
fclose(fileID);
end

end