function surrogate_test(y, surrogate_algorithm, Nsurr, permutation_order, out_name, param)
%%%%%%%%%%%%   Test if data are deterministic  %%%%%%%%%%%%%%%%%%%%%
disp(' Test if data are deterministic');
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
%out_name = append('output_surrogate_analysis_piecewise.txt'); 
fileID = fopen(out_name,'a');
fprintf(fileID, '%d, %d, %d, %d, %d, %d, %d, %d \n',param, perm_h, ...
        min(surr_h), max(surr_h), length(y), permutation_order, Nsurr, step2_time);