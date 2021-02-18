clear all
surrogate_algorithm = 'AAFT'
Nsurr = 100
permutation_order = 5
no_data = 2e5
% load data
for mV = 0.1%[-200 -160:20:80 200]
if mV==200
    name = append('../../SUBSTATES_200/+200.txt');
elseif mV==-200
    name = append('../../SUBSTATES_200/-200.txt');
elseif mV == 0
    name = append('../../SUBSTATES/zero.txt');
elseif mV < 0
    name = append('../../SUBSTATES/',int2str(abs(mV)) ,'.txt');
elseif mV > 0
    name = append('../../SUBSTATES_positive/',int2str(abs(mV)) ,'.txt');
end
%name = '../../Markov_Models/3state_model/sample_states.txt';
%name = '../../Markov_Models/3state_model/signal_noise1.txt';
%name = '../../Markov_Models/3state_model/signal_noise2.txt';
%name = '../../Markov_Models/3state_model_indistingushable/signal_noise1.txt'
%name = '../../Markov_Models/3state_model_indistingushable/signal_noise2.txt'
 
%name = append('../../SUBSTATES_200/-60.txt')
%name = append('../../SUBSTATES_200/+60.txt')
%name = '../../SUBSTATES_200/-60mV_lowPo.txt' %It is incorporated the other ...
        %way round in the bilayer so that you need to invert it.
%name = append('../../SUBSTATES_200/+100.txt')
%name = '../../SUBSTATES_200/+60mV-NFX.txt'

%name = '../../BetaFunctionAnalysis/TREK-2_4-state_model_STAR-SHAPED_fit/4-state_starModel_downsampled_noNoise.txt'
%name = '../../BetaFunctionAnalysis/TREK-2_4-state_model_STAR-SHAPED_fit/4-state_starModel_downsampled_NoiseAfterFilter.txt'
%name = '../../BetaFunctionAnalysis/TREK-2_4-state_model_STAR-SHAPED_fit/4-state_starModel_downsampled_NoiseBeforeFilter.txt'
%name = '../../BetaFunctionAnalysis/FITS/-200mV/-200_mV_noNoise.txt'
name = '../../BetaFunctionAnalysis/FITS/-200mV/-200_mV_Noise.txt'
table = readtable(name);
conv =  table2array(table);
disp('size of array (converted from table)');
disp(size(conv));
clear table;
 
subset = min(no_data, size(conv,1))
if mV == 0.1 
    y = conv(1:subset,1);
elseif mV ==0.2
    y = conv(1:subset,1);
else
    y = -1*conv(1:subset,2);
end
clear conv

% Step 2: Test if data are deterministic
    disp('Step 2: Test if data are deterministic');
    step2 = tic();
    surr_y=y;

    if size(surr_y,1)>size(surr_y,2)
    	disp('transpose vector')
        surr_y=surr_y';
    end
    
    % Test if data are stochastic and linear using permutation entropy
    % and AAFT surrogates
    
        try
            [surr, params] = surrogate(surr_y, Nsurr, surrogate_algorithm, 1, 1);
        catch
            disp('surroAAgate analysis with zscore!')
            [surr, params] = surrogate(zscore(surr_y), Nsurr, surrogate_algorithm, 1, 1);
        end
        sig=params.cutsig;
        for i = 1:Nsurr
            surr_h(i) = petropy(surr(i,:),permutation_order,1);
        end
        perm_h = petropy(sig,permutation_order,1);
        
        stochastic = perm_h>=min(surr_h)&&perm_h<=max(surr_h)
        disp('Permutation entropy trajectory:')
        disp(perm_h)
        disp('Min/max of surrogates')
        disp([min(surr_h)  max(surr_h) ])
  
    step2_time = toc(step2)
   %out_name = append('output_surrogate_analysis.txt'); 
   %out_name = append('output_surrogate_analysis_NFX.txt'); 
   %out_name = append('output_Markov_surrogate_analysis.txt');
   out_name = append('output_MarkovBessel_surrogate_analysis.txt');
    fileID = fopen(out_name,'a');
    fprintf(fileID, '%d, %d, %d, %d, %d, %d, %d, %d \n',mV, perm_h, ...
        min(surr_h), max(surr_h), subset, permutation_order, Nsurr, step2_time);
end