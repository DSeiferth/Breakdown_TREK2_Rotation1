clear all, clc, close all
%mV, perm_h, min(surr_h), max(surr_h), subset, permutation_order, Nsurr, step2_time

out_name = 'output_surrogate_analysis.txt';
table = readtable(out_name);
mV = table2array(table(:,1));
H = table2array(table(:,2));
surr_min = table2array(table(:,3));
surr_max = table2array(table(:,4));
subset = table2array(table(:,5));
perm_order = table2array(table(:,6));
Nsurr = table2array(table(:,7));

order = 5
max_H = log2(factorial(order))

% low and high open prob
out_name = 'output_surrogate_analysis_low_high_Po.txt';
table = readtable(out_name);
H_lowPo_minus60 = table2array(table(3,2))/max_H; 
H_highPo_plus60 = table2array(table(2,2))/max_H; 
H_highPo_minus60 = table2array(table(1,2))/max_H; 
H_highPo_plus100 = table2array(table(4,2))/max_H; 


figure(1);
plot(mV, H/max_H,'-x', 'color', 'b')
hold on
plot(0, H(10)/max_H, '-s','MarkerSize',10,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6])
hold on
plot(mV, surr_min/max_H,'-x', 'color', 'red')
hold on
plot([-60], H_lowPo_minus60, '-s','MarkerSize',7,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor','black')
hold on 
plot([-60], H_highPo_minus60, '-s','MarkerSize',7,...
    'MarkerEdgeColor','green',...
    'MarkerFaceColor','green')
hold on
plot([60], H_highPo_plus60, '-s','MarkerSize',7,...
    'MarkerEdgeColor','green',...
    'MarkerFaceColor','green')
hold on
plot([100], H_highPo_plus100, '-s','MarkerSize',7,...
    'MarkerEdgeColor','green',...
    'MarkerFaceColor','green')
legend('H(original data) -> chaotic', 'H(membrane noise) -> stochastic','min H(surrogates)', ...
    'H(low open prob.) -> chaotic', 'H(high open prob.) -> chaotic' ,...
    'Location', 'NorthEast')
xlabel('Membrane potential [mV]')
ylabel('Permutation entropy H(n)/log(n!)')
%plotname = 'plot_permutation_entropy';
plotname = 'plot_permutation_entropy_surrogates_included';
title({'Permutation order n = 5','sampling step = 5\mus', 'AAFT surrogates'})
set(gca,'FontSize',14)
%set(gca, 'XScale', 'log')
saveas(gcf,plotname,'epsc');