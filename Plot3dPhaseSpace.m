function Plot3dPhaseSpace(x_t,ndim, tau, subset, N_scatter, plotname, title_name, title_name2)
datcnt = subset;
ires = 10;
maxbox = 5e4 %50000 until 3e5; %7664
dummy = 'dummy.dat';
writematrix(x_t, dummy); 
db = basgen(dummy, tau, ndim, ires, datcnt, maxbox);
dataplot = [];
delay = 0:tau:(ndim-1)*tau;
data = db.data;
for ii = 1:(datcnt-(ndim-1)*tau)
        dataplot = [dataplot; data(ii+delay)];
end

fh1 = figure(1);
ax1 = axes(fh1);
grid(ax1, 'on');
set(fh1, 'Renderer', 'OpenGL');
%plot3(dataplot(:,1), dataplot(:,2), dataplot(:,3), '.', 'MarkerSize', 3)
plot3(dataplot(1:N_scatter,1), dataplot(1:N_scatter,2), dataplot(1:N_scatter,3), '.', 'MarkerSize', 3)
title({ 
    [title_name]
    [title_name2]
    [' embedding: m = ' num2str(ndim) '; \tau = ' num2str(tau)]
});
xlabel('x(i)')
ylabel('x(i+\tau)')
zlabel('x(i+2\tau)')
grid(ax1, 'on');
set(gca,'FontSize',15)
saveas(gcf,plotname,'epsc');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% https://uk.mathworks.com/matlabcentral/fileexchange/23897-n-dimensional-histogram
% Like histc but for n-dimension
Nk = 500
[count edges mid loc] = histcn(dataplot, Nk);
% make a grid for plotting
[X Y Z]=ndgrid(edges{1}, edges{2}, edges{3});
X=X(:); Y=Y(:); Z=Z(:);
% logscale for size
log_count = log(count(:));
log_count(log_count == -Inf)=0;

% calculate sizes so the most dense cell gets a value of 100
% also convert from volume to "area" (as if drawing a sphere with
% the right volume and cross-sectional area s)
maxVal = 10
s_scale = maxVal/(max(log_count(:))^(2/3));
s = log_count(:).^(2/3) * s_scale;
% convert any zeros to small numbers for scatter3
s(s==0)=realmin;
%https://uk.mathworks.com/matlabcentral/answers/502563-how-to-plot-a-scattered-heat-map
%color map
cmap = jet(256);
bla = log(count(:));
bla(bla == -Inf)=0;
v = rescale(bla, 1, 256);
numValues = length(count(:))
markerColors = zeros(numValues, 3);
% Now assign marker colors according to the value of the data.
for k = 1 : numValues
    row = round(v(k));
    markerColors(k, :) = cmap(row, :);
end
% plot the densities
fh = figure(2);
set(fh, 'Renderer', 'OpenGL'); % faster drawing
scatter3(X, Y, Z, s, markerColors, 'filled'); 
title({ 
    [title_name]
    [title_name2]
    [' embedding: m = ' num2str(ndim) '; \tau = ' num2str(tau)]
});
xlabel('x(i)')
ylabel('x(i+\tau)')
zlabel('x(i+2\tau)')
set(gca,'FontSize',15)
saveas(gcf,append(plotname, '_hist'),'epsc');
end

