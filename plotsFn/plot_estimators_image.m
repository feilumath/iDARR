function plot_estimators_image(est_array,tgrid,L_operator,y,y_true,lgnd,new_figure)
% plot the image of the estimators

colors; % get linestyle and dark colors

if ~exist('new_figure','var'); new_figure=1; end
if new_figure==1; figure; end 
[n_vec,n_est] = size(est_array); 

plot(tgrid, y,':x','linewidth',0.5,'Color',[0.7,0.85,0.50]);hold on;
plot(tgrid, y_true,'-k','linewidth',1);hold on;
for n=1:n_est-1
    f_est = est_array(:,n);
    plot(tgrid, L_operator*f_est,linestyle{n+1},'Color',dred_do_db(n,:),'linewidth',1);hold on;
end
hold off;
xlim([tgrid(1) tgrid(end)])
lgnd = ['Observed',lgnd(1:end),];
title('Estimated y');legend(lgnd);
xlabel('t'); ylabel('y(t)');
end