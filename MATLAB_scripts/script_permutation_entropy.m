for mV = [-200, 200]%-160:20:80
disp(append('Voltage = ',int2str(mV)));
if mV==200
    name = append('../../SUBSTATES_200/+200.txt');
elseif mV==-200
    name = append('../../SUBSTATES_200/-200.txt');
elseif mV<0
    name = append('../../SUBSTATES/',int2str(-mV) ,'.txt');
elseif mV == 0
    name = append('../../SUBSTATES/zero.txt');
elseif mV > 0
    name = append('../../SUBSTATES_positive/',int2str(mV) ,'.txt');
end
disp(name);
table = readtable(name);
conv =  table2array(table);
disp('size of array (converted from table)');
disp(size(conv));
clear table;

for subset =  2.8e5%[logspace(2,6), size(conv,1)]
t = conv(1:subset,1); % column
X = conv(1:subset,2);
disp('size of X');
disp(size(X));

% H = pentropy(x,n,tau,method,accu)
for n = 5; % permutation order
for tau = 1%[1, 5, 10, 20, 100]; % time lag 
now1 = tic;
H = pentropy(X,n,tau);
disp(H);
wholeTime = toc(now1)

% write output in file
%out_name = append('output_permutation_entropy_order',int2str(n),'.txt');
%out_name = append('output_permutation_entropy_timelag',int2str(tau),'.txt');
out_name = 'output_permutation_entropy_200.txt'
fileID = fopen(out_name,'a');
fprintf(fileID, '%d, %d, %d, %d, %d, %d \n',mV, H, n,subset, wholeTime, tau);
fclose(fileID);

end
end
end
end