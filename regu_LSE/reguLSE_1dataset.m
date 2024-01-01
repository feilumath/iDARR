%% computing estimator and its errors 
function [err_L2rho,errors,err_l2,loss_array,iter_stop,dim_Krylov] = reguLSE_1dataset(A,b,Abar,bbar,B,f_true,V_AB,V_A,r,xgrid,rho,eigAB,plotON,normType,method)
    if ~exist('plotON','var'); plotON=1;end
    % plotON = 1; 
    %% test different regularizations: l2, L2, RKHS, and may be H1? 
    switch method
        case 'Direct'  % DARTR 
            [~,~,cLcurveAll] = reg_Lcurve_3in1(A,b,B,plotON,normType); 
            est_array        = [cLcurveAll.creg_l2,cLcurveAll.creg_L2,cLcurveAll.creg_RKHS,f_true];
            iter_stop        = [cLcurveAll.creg_l2_lambda_opt,cLcurveAll.creg_L2_lambda_opt,cLcurveAll.creg_RKHS_lambda_opt];
            dim_Krylov       = rank(A)*ones(1,3);
            [errors]         = compute_estimator_Error(B,V_AB,eigAB,V_A,r,est_array,rho,xgrid); 
            err_L2rho        = errors.L2rho_l2_projAB(1,:);
            err_l2           = errors.L2rho_l2_projA(2,:);
            loss_array       = [cLcurveAll.creg_l2_loss_val,cLcurveAll.creg_L2_loss_val,cLcurveAll.creg_RKHS_loss_val];

        case 'IR'
            [~,cLcurveAll] = IRreg_3in1(A,b,B,rho,plotON,normType); 
            iter_stop        = [cLcurveAll.creg_lsqr_k,cLcurveAll.creg_wlsqr_k,cLcurveAll.creg_dartr_k];
            dim_Krylov       = [cLcurveAll.creg_lsqr_k_terminate,cLcurveAll.creg_wlsqr_k_terminate,cLcurveAll.creg_dartr_k_terminate];
            est_array        = [cLcurveAll.creg_lsqr,cLcurveAll.creg_wlsqr,cLcurveAll.creg_dartr,f_true];
            [errors]         = compute_estimator_Error(B,V_AB,eigAB,V_A,r,est_array,rho,xgrid); 
            err_L2rho        = errors.L2rho_l2_projAB(1,:);
            err_l2           = errors.L2rho_l2_projA(2,:);
            loss_array       = [cLcurveAll.creg_lsqr_loss_val,cLcurveAll.creg_wlsqr_loss_val,cLcurveAll.creg_dartr_loss_val];
    end
end