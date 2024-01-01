function [err_L2rho_projAB,err_cells,err_l2_projA,loss_array,iter_stop] = multi_simuls(f_true,Abar,B,V_AB,V_A,r,xgrid,rho,eigAB,sysInfo,nsr_seq,normType,method,saveDIR)
% if exist(saveDIR,'file')
%     load(saveDIR); return; 
% end
% multiple simulations: 
n_simuls = 100; 
% normType  = {'l2','L2','RKHS'}; 
err_L2rho_projAB = zeros(n_simuls,length(normType),length(nsr_seq)); 
err_cells = cell(n_simuls,length(nsr_seq)); 
plotON    =0; 

A = sysInfo.L_operator;
%% nsr = 0: only 1 dataset since there is no noise;
nsr = 0; 
tic
b_true = A*f_true; 
b_true_norm = sqrt(sum(b_true.^2)*sysInfo.dt/sysInfo.T);
% b      = b_true + b_true_norm*nsr*randn(sysInfo.tn, 1)*sqrt(sysInfo.dt);
b      = b_true + b_true_norm*nsr*randn(sysInfo.tn, 1);
bbar      = A'*b;
[err_L2rho_projAB_nrs0,errors_nrs0] = reguLSE_1dataset(A,b,Abar,bbar,B,f_true,V_AB,V_A,r,xgrid,rho,eigAB,plotON,normType,method);
t1  =toc; 
fprintf('nsr = 0 with 1 dataset. Time elapsed: %2.2f minutes \n', t1/60); 
fprintf('Total time to expect: %2.2f minutes \n ', t1/60*n_simuls* length(nsr_seq));


err_l2_projA   = err_L2rho_projAB; 
loss_array     = err_L2rho_projAB; 
iter_stop      = err_L2rho_projAB;
dim_Krylov     = err_L2rho_projAB; 
%% nsr > 0: 100 datasets   
for m = 1:length(nsr_seq) 
    nsr = nsr_seq(m);  tn =sysInfo.tn;
    parfor n=1:n_simuls
        b      = b_true + b_true_norm*nsr*randn(tn, 1);
        bbar      = A'*b;
        [err_L2rho_projAB(n,:,m), err_cells{n,m},err_l2_projA(n,:,m),loss_val,k_iter,dimK] = ...
            reguLSE_1dataset(A,b,Abar,bbar,B,f_true,V_AB,V_A,r,xgrid,rho,eigAB,plotON,normType,method);
         loss_array(n,:,m)       = loss_val; 
         iter_stop(n,:,m)        = k_iter; 
         dim_Krylov(n,:,m)       = dimK; 
    end
    fprintf('Progress %i out of %i \n ',  m,length(nsr_seq)); 
end
fprintf('DONE \n'); 

if exist('saveDIR','var') 
    save(saveDIR,'err_cells','err_L2rho_projAB','err_L2rho_projAB_nrs0',...
        'errors_nrs0','err_l2_projA','loss_array','iter_stop','dim_Krylov'); 
end

end

