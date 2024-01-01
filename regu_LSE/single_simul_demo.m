%% Get regression vector for f. 
switch kernel_type
    case 'exp'
        nsr = 0.5;          % noise to signal ratio
    case 'poly'
        nsr   = 0.0625;    % noise to signal ratio    
end           

if ~exist('f_true','var')
   f_true_func = @(x) x.^2;
   f_true = f_true_func(xgrid);
end

f_true = f_true /(dx*sum(f_true));

sqrt_deltat = sqrt(sysInfo.dt); 
b_true = A*f_true; 
b_true_norm = sqrt(sum(b_true.^2)*sysInfo.dt/sysInfo.T);
b      = b_true + b_true_norm*nsr*randn(tn, 1);

bbar        = A'*b;       % Abar, bbar are used in DARTR.  
bbar_true   = A'*b_true;

%% unconstained LSE with regularuzations: l2, L2, RKHS
plotON = 1;  
if exist('fig_dir','var') 
    figname = [fig_dir,file_str,'_estimator']; 
else 
    clear figname; 
end
%% test different regularizations: l2, L2, RKHS, and may be H1? 
DnormType  = {'l2','L2','RKHS'}; 
IRnormType  = {'IR-l2','IR-L2','iDARR'};

switch method
    case 'Direct'
         [~,~,cLcurveAll] = reg_Lcurve_3in1(A,b,B,plotON,DnormType);  % L-curve method
    case 'IR'
        [~,cLcurveAll] = IRreg_3in1(A,b,B,rho,plotON,IRnormType);  % iterative method
end


 
%%
if plotON==1
    switch method
        case 'Direct'
            lgnd      = {'True','Est-l2','Est-L2','Est-RKHS'};
            est_array = [cLcurveAll.creg_l2,cLcurveAll.creg_L2,cLcurveAll.creg_RKHS,f_true];
 
        case 'IR'
            lgnd      = {'True','IR-l2','IR-L2','iDARR'};
            est_array = [cLcurveAll.creg_lsqr,cLcurveAll.creg_wlsqr,cLcurveAll.creg_dartr,f_true];       
    end
    new_figure = 0;
    h= figure;
    subplot(121); plot_estimators(est_array,lgnd,xgrid,rho,new_figure);
    subplot(122); plot_estimators_image(est_array,tgrid,A,b,b_true,lgnd,new_figure);   
    if exist('figname','var')
       figure(h); set_positionFontsAll;   % print([figname,'estimator.pdf'],'-dpdf', '-bestfit');
    end
    
    figure; 
    subplot(121);
    % projection to the V_AB
    titl = 'EigenAB projection';
    plot_estimators_projection(B,V_AB,min(2*r,length(B(1,:))),est_array,lgnd,titl,new_figure);
    % project to V_A
    subplot(122);
    titl = 'EigenA projection';
    plot_estimators_projection(eye(xn),V_A,min(2*r,length(B(1,:))),est_array,lgnd,titl,new_figure);

    switch method
        case 'IR'
            h= figure;
            res_array = {cLcurveAll.creg_lsqr_res,cLcurveAll.creg_wlsqr_res,cLcurveAll.creg_dartr_res};
            k_array   = [cLcurveAll.creg_lsqr_k,cLcurveAll.creg_wlsqr_k,cLcurveAll.creg_dartr_k];
            fig_name = [figname,'_res'];
            plot_res(res_array,k_array,lgnd,fig_name,new_figure);
    end
end

[errors] = compute_estimator_Error(B,V_AB,eigAB,V_A,r,est_array,rho,xgrid); 

display(errors.L2rho_l2_discription)
display(errors.L2rho_l2)


display(errors.discription_L2rho_l2_projAB)
display(errors.L2rho_l2_projAB)

display(errors.discription_L2rho_l2_projA)
display(errors.L2rho_l2_projA)

switch method
    case 'Direct'
        fprintf('Loss values: l2,L2,RKHS \n')
        loss_array  = [cLcurveAll.creg_l2_loss_val,cLcurveAll.creg_L2_loss_val,cLcurveAll.creg_RKHS_loss_val]

    case 'IR'
        fprintf('Loss values: IR l2,IR L2,iDARR \n')
        loss_array  = [cLcurveAll.creg_lsqr_loss_val,cLcurveAll.creg_wlsqr_loss_val,cLcurveAll.creg_dartr_loss_val]
        fprintf('Dimension of Krylov subspace and k_optimal: IR l2,IR L2,iDARR \n')
        dim_Krylov_subspace = [cLcurveAll.creg_lsqr_k_terminate,cLcurveAll.creg_wlsqr_k_terminate,cLcurveAll.creg_dartr_k_terminate; 
            cLcurveAll.creg_lsqr_k,cLcurveAll.creg_wlsqr_k,cLcurveAll.creg_dartr_k ]
end
