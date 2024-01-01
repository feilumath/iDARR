function h = plot_mean_std(x_ind,data_array,newfigure,lgnd,label_y,label_x,option,figname) 
% plot the ensemble mean with std bar of samples in data_array:  X(samples,types,x_ind)
%{
    x_ind      - the index variable for the data 
    data_array -  n_samples x n_types x n_ind
%}
[n_samples,n_types,n_ind] = size(data_array); 
std_data    = zeros(n_types,n_ind);
mean_data   = zeros(n_types,n_ind);
for nn= 1:n_ind
    temp     = squeeze(data_array(:,:,nn));
    std_data(:,nn) = std(temp,0,1);
    mean_data(:,nn)= mean(temp,1);
    % boxplot(temp,'Notch','on');
end
fprintf('\n Mean of %i simulations \n ',n_samples);
disp(mean_data);


if exist('newfigure','var') && newfigure ==1; figure; end 
colors;     % colors and linestyples 

switch option
    case 'log'
        mean_data   = log10(mean_data); 
        data_array  = log10(data_array);
end

for i= 1:n_types
    h(i) = boxchart(squeeze(data_array(:,i,:)),'BoxWidth',0.3,'BoxFaceColor',dred_do_db(i,:), ...
        'BoxEdgeColor',dred_do_db(i,:),'WhiskerLineColor',dred_do_db(i,:),'JitterOutliers','on',...
        'MarkerStyle',markerstyle(i),'MarkerColor',dred_do_db(i,:),'MarkerSize',3); hold on;
    p(i) = plot(mean_data(i,:),linestyles{i},'linewidth',2,'Color',dred_do_db(i,:),'MarkerSize',10); hold on;
end
hold off
xticklabels(string(x_ind));


box on
xlabel(label_x); 

switch option
    case 'log'
        label_y = ['log_{10} ', label_y];
end
if exist('lgnd','var'); legend(p,lgnd,'Location','best');  end
if exist('label_y','var'); ylabel(label_y); yticklabels('auto'); end
if exist('figname','var')
    set_positionFontsAll;   
     print([figname,'.pdf'],'-dpdf','-bestfit'); 
end
end


