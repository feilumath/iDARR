%% 
nsr_seq  = [0.0625,0.125,0.25,0.5,1];     % noise to signal ratio    
DnormType  = {'l2','L2','RKHS'}; 
IRnormType  = {'IR-l2','IR-L2','iDARR'};


%% Part I: f_true outside the x) (sin(x-6)).^2+1;   % True f
f_true_func = @(x) x.^2;
f_true = f_true_func(xgrid);


figure; plot(xgrid,f_true);  %/(dx*sum(f_true))); 

%
file_str  =['outsideFSOI_',method,kernel_type]; 
data_name = [SAVE_DIR,'/data_',file_str,'.mat']; 
fig_dir  = [SAVE_DIR,'/figures/']; if ~exist(fig_dir,'dir'), mkdir(fig_dir); end  

%% Part I.1  single simulultion demo and tuning. 
single_simul_demo; 

%% Part I.2: multiple simulations 
fprintf('Multiple tests: f true outside the FSOI\n ');

switch method
    case 'Direct'
        normType = DnormType;
    case 'IR'
        normType = IRnormType;  % iterative method
end


[err_L2rho_projAB,err_cells,err_l2_projA,loss_array,iter_stop] = multi_simuls(f_true,Abar,B,V_AB,V_A,r,xgrid,rho,eigAB,sysInfo,nsr_seq,normType,method,data_name);


newfigure =1; 
plot_scale = 'log';  % 'log' or 'linear'

label_y = 'L^2(\rho) error';   label_x = 'nsr'; 
fig_name = [fig_dir,file_str,'_L2rho_error']; 
plot_mean_std(nsr_seq,err_L2rho_projAB,newfigure,normType,label_y,label_x,plot_scale,fig_name);  

label_y = 'l2 error';
fig_name = [fig_dir,file_str,'_l2error']; 
plot_mean_std(nsr_seq,err_l2_projA,newfigure,normType,label_y,label_x,plot_scale,fig_name);  

label_y = 'Loss value';
fig_name = [fig_dir,file_str,'_loss']; 
plot_mean_std(nsr_seq,loss_array,newfigure,normType,label_y,label_x,plot_scale,fig_name);  

plot_scale = 'linear';
label_y = 'Stop Iteration';
fig_name = [fig_dir,file_str,'_stopiter']; 
plot_mean_std(nsr_seq,iter_stop,newfigure,normType,label_y,label_x,plot_scale,fig_name);  


%%  Part II: f_true inside the FSOI: 
f_true = V_AB(:,2);   % f_true = V_AB(:,5);    % ind>5 ( i.e, eigAB<1e-9): rkhs not good, l2+L2 can slightly tolerate more, bc. not using rkhs inversion).
f_true = f_true/(dx*sum(f_true)); 

%%
file_str  = ['insideFSOI_',method,kernel_type]; 
data_name = [SAVE_DIR,'/data_',file_str,'.mat']; 
fig_dir  = [SAVE_DIR,'/figures/']; if ~exist(fig_dir,'dir'), mkdir(fig_dir); end  

%% Part II.1: single simulultion demo 
single_simul_demo; 

%%  Part II.2: multiple simulations 
fprintf('Multiple tests: f true inside the FSOI\n ');

switch method
    case 'Direct'
        normType = DnormType;
    case 'IR'
        normType = IRnormType;  % iterative method
end

[err_L2rho_projAB2,err_cells,err_l2_projA2,loss_array,iter_stop] = multi_simuls(f_true,Abar,B,V_AB,V_A,r,xgrid,rho,eigAB,sysInfo,nsr_seq,normType,method,data_name);

newfigure =1;
plot_scale = 'log';  % 'log' or 'linear'

label_y = 'L^2(\rho) error';   label_x = 'nsr'; 
fig_name = [fig_dir,file_str,'_L2rho_error']; 
plot_mean_std(nsr_seq,err_L2rho_projAB2,newfigure,normType,label_y,label_x,plot_scale,fig_name); 

label_y = 'l2 error';
fig_name = [fig_dir,file_str,'_l2error']; 
plot_mean_std(nsr_seq,err_l2_projA2,newfigure,normType,label_y,label_x,plot_scale,fig_name);  

label_y = 'Loss value';
fig_name = [fig_dir,file_str,'_loss']; 
plot_mean_std(nsr_seq,loss_array,newfigure,normType,label_y,label_x,plot_scale,fig_name);  

plot_scale = 'linear';
label_y = 'Stop Iteration';
fig_name = [fig_dir,file_str,'_stopiter']; 
plot_mean_std(nsr_seq,10.^(iter_stop),newfigure,normType,label_y,label_x,plot_scale,fig_name);  

