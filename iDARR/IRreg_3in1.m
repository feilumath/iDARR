function [x_reg,cLcurveAll] = IRreg_3in1(A,b,~,rho,plotON, normType)
% regularizaiton using iterative methods, considering three norms: l2, L2, RKHS
% NOTE: IRhybrid_lsqr Hybrid version of LSQR algorithm
%{
normType  = {'l2','L2','RKHS'} or part of them; other types are possible. 
                  x_l = argmin_x  |Ax -b|^2 + lambda x'*B*x
                      = (A'*A+lambda*B)\(A'*b)  -- not using this LSE! 
use RegParam to set  a method to find the regularization, 'gcv', 'wgcv' or 'modgcv'
Input:  
        - A, b: Ax = b;   
        - rho: exploration measure, = sum(A,1) if not specified. 
        - B: basis matrix.  2023-10-17: IR-DA-RKHS works for basis matrix B = diag(\rho) 
Regulariztaion methods:
    norm_type = 'l2':     |Ax -b|^2 + lambda x'* Id *x
    norm_type = 'L2':     |Ax -b|^2 + lambda x'* B *x,   with B being the L2 norm matrix
    norm_type = 'RKHS':   |Ax -b|^2 + lambda x'* B-rkhs *x,
% Copyright(c): Fei Lu, Jinchao Feng, Haibo Li 
%}

% Avoid A2, b2: They are used in direct matrix-decomposition based method.  
% A2 = sqrtm(Abar);   b2 = pinv(A2)*bbar; % Costly! Avoid. 

%% TBD: 
% Issue1: how do we know k=10? 

% B = 0*B; % not used in this version yet. TBD later.  
            
max_iter = length(A(1,:))-1; 


%% Loop for all normTypes 
for nn = 1:length(normType)
    norm_type = normType{nn};
    switch norm_type
        case 'IR-l2'
            k = max_iter;
            reorth = 1;
            [X1, res1, eta1] = lsqr_b(A, b, k, reorth);
            [k1_cor,info1] = Lcurve(res1,eta1,plotON,'l2-l2');

            x_reg                            = X1(:,k1_cor);
            cLcurveAll.creg_lsqr             = x_reg;
            % cLcurveAll.creg_lsqr_loss_val    = norm(A2*x_reg - b2);
            cLcurveAll.creg_lsqr_loss_val    = norm(A*x_reg - b);  % l2 residual norm^2 = the loss value
            cLcurveAll.creg_lsqr_k           = k1_cor;
            k_terminate      = length(eta1);
            cLcurveAll.creg_lsqr_k_terminate = k_terminate; 
            cLcurveAll.creg_lsqr_res     = res1;

        case 'IR-L2'
            k = max_iter;
            reorth = 1;
            [X2, res2, eta2] = wlsqr(A, rho, b, k, reorth);
            [k2_cor,info2]   = Lcurve(res2,eta2,plotON,'l2-L2');
          
            x_reg                             = X2(:,k2_cor);
            cLcurveAll.creg_wlsqr             = x_reg;
            cLcurveAll.creg_wlsqr_loss_val    = norm(A*x_reg - b);
            cLcurveAll.creg_wlsqr_k             = k2_cor;

             k_terminate      = length(eta2);
            cLcurveAll.creg_wlsqr_k_terminate = k_terminate;
            cLcurveAll.creg_wlsqr_res     = res2;

        case 'iDARR'
            k = max_iter;
            reorth = 1;
            [X3, res3, eta3] = dartr_spr(A, b, rho, k, reorth);
            [k3_cor,info3]   = Lcurve(res3,eta3,plotON,'l2-xHG');
%             if isempty(k3_cor)
%                 k3_cor     =  find(log10(res3)>-10,1,'last') -1; 
%             end
            x_reg                             = X3(:,k3_cor);
            cLcurveAll.creg_dartr             = x_reg;
            cLcurveAll.creg_dartr_loss_val    = norm(A*x_reg - b);
            cLcurveAll.creg_dartr_k           = k3_cor;
            k_terminate      = length(eta3);
            cLcurveAll.creg_dartr_k_terminate = k_terminate; 
            cLcurveAll.creg_dartr_res     = res3;
    end
end

end



