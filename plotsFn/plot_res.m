function plot_res(res_array,k_array,lgnd,figname,new_figure)
% plot the residuals 
% Input 
%     Est_array: each cell is a vector of residuals


colors; % get linestyle and dark colors

n_cell = length(res_array); 
kterminate = zeros(n_cell,1);
for n = 1:n_cell
    kterminate(n) = length(res_array{n});
end
if ~exist('new_figure','var'); new_figure=1; end
if new_figure==1; figure; end 

for n=1:n_cell
    p(n) = plot(1:kterminate(n),log10(res_array{n}),linestyle{n+1},'Color',dred_do_db(n,:),'linewidth',1.5); hold on;
    plot(k_array(n),log10(res_array{n}(k_array(n))),'Marker',markerstyle(n),'MarkerEdgeColor',dred_do_db(n,:),'MarkerSize',12,'linewidth',2); hold on;
end
 
xlabel('iteration'); ylabel('log_{10} ||Ax_k-b||_2'); 
lgnd = lgnd(2:end); 
legend(p,lgnd,'Location','best');

if exist('figname','var')
    set_positionFontsAll;   
    print([figname,'.pdf'],'-dpdf','-bestfit'); 
end

end
