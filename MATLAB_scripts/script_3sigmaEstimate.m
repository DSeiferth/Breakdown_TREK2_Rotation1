clear all, close all, clc
cutOffFreq = 10
folder = '5kHz/';
step = 5e-6;
for mV = -200
save=1
fixedEdges = -80:0.33:10;

subset = 1e4;
disp(append('Voltage = ',int2str(mV)));
name = append('../SUBSTATES_differentFreq/', num2str(cutOffFreq), 'kHz', num2str(mV),'mV.txt');

table = readtable(name);
conv =  table2array(table);
close all
clear table;
x_1 = conv(:,2);

[N,edges] = histcounts(x_1, fixedEdges,  'Normalization','pdf');
edges = edges(2:end) - (edges(2)-edges(1))/2;
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
title({[int2str(mV), 'mV ',  'cut-off freq. ', num2str(cutOffFreq), 'kHz']
    ['channel in closed state - noise estimation']
    }, 'FontWeight','Normal')
xlabel('Time [s]')
ylabel('Re-scaled Current')
set(gca,'FontSize',14)
plotname = append(folder, 'noiseLevel_',int2str(mV) ,'mV_', num2str(cutOffFreq), 'kHz');
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
title({
    [ 'cut-off freq. ', num2str(cutOffFreq), 'kHz']
    ['Noise of closed state, \sigma = ', num2str( stdClosed )]
    }, 'FontWeight','Normal')
set(gca,'FontSize',14)
plotname = append(folder, 'noiseLevel_distr_',int2str(mV) ,'mV_', num2str(cutOffFreq), 'kHz');
if save == 1
saveas(gcf,plotname,'epsc');
end
%%%%%%%%%%%%   Output  modes and sigma  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save == 1
out_name = append(folder, 'mV_cutOffFreq_Peak1_Peak2_stdClosed_Po.txt');
fileID = fopen(out_name,'a');
fprintf(fileID, '%d, %d, %d, %d, %d, %d\n',mV, cutOffFreq, posMode(1),posMode(2),stdClosed, open);
fclose(fileID);
end

%%%%%%%%%%%%%  Threshold criterium   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
threshold = posMode(2)-3*stdClosed;
indxClosed = x_1 > threshold;
indxOpen = x_1 < threshold;
IdealisedTS = x_1;
IdealisedTS(indxClosed) = posMode(2); % closed level
IdealisedTS(indxOpen) = posMode(1); % open level

dur = 0;
prev = x_1(1);
duration = [];
OpenClosed = [];
for i = 2:1:length(x_1)
    if (prev > threshold && x_1(i) > threshold) || (prev < threshold && x_1(i) < threshold)
        dur = dur + step;
    else
        duration(end+1) = dur;
        dur = 0;
        if prev > threshold % closed
            OpenClosed(end+1) = 0;
        elseif prev < threshold % open
            OpenClosed(end+1) = -1;
        end
        prev = x_1(i);
    end      
end

if save == 1
    fname1 = append(folder,  num2str(cutOffFreq), 'kHz_', num2str(mV),'mV_OpenClosed.txt');
    writematrix(OpenClosed', fname1);
    fname1 = append(folder,  num2str(cutOffFreq), 'kHz_', num2str(mV),'mV_duration.txt');
    writematrix(duration', fname1);
end
%%%%%%%%%%%%   Plot dwell time histogram  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5)
idxOpen = OpenClosed == -1;
idxClosed = OpenClosed == 0;
[N_open,edges_open] = histcounts(duration(idxOpen), 'Normalization','pdf');
edges_open = edges_open(2:end) - (edges_open(2)-edges_open(1))/2;
plot(edges_open, N_open);
hold on
[N_closed,edges_closed] = histcounts(duration(idxClosed), 'Normalization','pdf');
edges_closed = edges_closed(2:end) - (edges_closed(2)-edges_closed(1))/2;
plot(edges_closed, N_closed);
legend('open', 'closed')
ylabel('pdf')
xlabel('Dwell time [s]')
%set(gca, 'YScale', 'log')
title({[int2str(mV), 'mV; ',  'cut-off freq. ', num2str(cutOffFreq), 'kHz']
    ['P_o = ', num2str(open), '; \sigma = ', num2str(stdClosed)]
    }, 'FontWeight','Normal')
set(gca,'FontSize',14)
plotname = append(folder, 'DwellTime_',int2str(mV) ,'mV_', num2str(cutOffFreq), 'kHz');
if save == 1
    saveas(gcf,plotname,'epsc');
end

figure(6)
log_edges = logspace(-6, -2, 25);
[N_open,edges_open] = histcounts(duration(idxOpen), log_edges, 'Normalization','count');
[N_closed,edges_closed] = histcounts(duration(idxClosed), log_edges, 'Normalization','count');
centered_log_edges = log_edges(1:end-1);
for i = 1:1:length(log_edges)-1
    centered_log_edges(i) = (log_edges(i)+log_edges(i+1))/2;
end
plot(centered_log_edges, sqrt(N_open), '-x');
hold on
plot(centered_log_edges, sqrt(N_closed), '-o');
legend('open', 'closed')
ylabel('\surd{counts}')
xlabel('Dwell time [s]')
set(gca, 'XScale', 'log')
title({[int2str(mV), 'mV; ',  'cut-off freq. ', num2str(cutOffFreq), 'kHz']
    ['P_o = ', num2str(open), '; \sigma = ', num2str(stdClosed)]
    }, 'FontWeight','Normal')
set(gca,'FontSize',14)
plotname = append(folder, 'DwellTime_log_',int2str(mV) ,'mV_', num2str(cutOffFreq), 'kHz');
if save == 1
    saveas(gcf,plotname,'epsc');
end

%%%%%%%%%%%%   Plot time series %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
step = 5e-6;
t_1 = [0:step:step*(subset-1)];
t_2 = [min(t_1), max(t_1)];
offset = 1e5;
plot(t_1, x_1(1+offset:subset+offset))
hold on
plot(t_2, [posMode(2), posMode(2)], 'k--')
hold on
plot(t_2, [posMode(2)-3*stdClosed, posMode(2)-3*stdClosed], 'k:', 'LineWidth', 2)
hold on
plot(t_1, IdealisedTS(1+offset:subset+offset), 'r')
xlabel('Time [s]')
ylabel('Current [pA]')
legend('Time series', 'Closed level', '3-\sigma boundary', 'Idealised time series')
title(append(int2str(mV), 'mV ', 'cut-off freq. ', num2str(cutOffFreq), 'kHz'), 'FontWeight','Normal')
set(gca,'FontSize',14)
plotname = append(folder, 'timeseries_',int2str(mV) ,'mV_', num2str(cutOffFreq), 'kHz');
if save == 1
saveas(gcf,plotname,'epsc');
end

%%%%%%%%%%%%%%%%    plot histograms        %%%%%%%%%%%%%%%%%%%%%%%%%%%
[N,edges] = histcounts(x_1, fixedEdges,  'Normalization','pdf');
edges = edges(2:end) - (edges(2)-edges(1))/2;
N_bounds = [0, max(N)];
figure(2)
plot(edges, N);
hold on 
plot([posMode(2), posMode(2)] , N_bounds, 'k--')
hold on
plot([posMode(2)-3*stdClosed, posMode(2)-3*stdClosed], N_bounds, 'k:', 'LineWidth', 2)
legend('Amplitude distribution', 'Closed level', '3-\sigma boundary', 'Location', 'NorthWest')
ylabel('pdf')
ylabel('Current [pA]')
%set(gca, 'YScale', 'log')
title({[int2str(mV), 'mV; ',  'cut-off freq. ', num2str(cutOffFreq), 'kHz']
    ['open mode at ', num2str(posOpenMode) ]
    ['closed mode at ', num2str(posClosedMode) ]
    ['P_o = ', num2str(open), '; \sigma = ', num2str(stdClosed)]
    }, 'FontWeight','Normal')
set(gca,'FontSize',14)
plotname = append(folder, 'distribution_',int2str(mV) ,'mV_', num2str(cutOffFreq), 'kHz');
if save == 1
saveas(gcf,plotname,'epsc');
fname1 = append(folder,  num2str(cutOffFreq), 'kHz_', num2str(mV),'mV_N.txt');
fname2 = append(folder,  num2str(cutOffFreq), 'kHz_', num2str(mV),'mV_edges.txt');
writematrix(N', fname1)
writematrix(edges', fname2)
end

end