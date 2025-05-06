%% iDARR demo: iterative Data-Adaptive RKHS regularization
% @c Haibo Li, Jinchao Feng, Fei Lu

add_mypaths_discrete;    % get SAVE_DIR = local dir for saving data


clc; close all; 

m= 500; n = 1000; 

expl_measure = 'L1';  % exploration measure
% expl_measure = 'L2';  % exploration measure

rng(2024)
% construct A: ill-conditioned 
A     = zeros(m,n);   %% We don't use random iid, because it is not ill-conditioned

xgrid = linspace(1, 5, n)';   dx = xgrid(2)- xgrid(1);
tgrid = linspace(0, 1, m)';
kernel_int = @(t, x) x.^(-1.).*abs(sin(x.*t +1)); % This phi is the integral kernel K(t,x)  ; 
                                                  % future: make it rough
% kernel_int = @(t, x) x.^(-1.).*abs(sin(x.*t +1))+ (x>2)+(x<4);                                                   
 for i = 1:m
    t = tgrid(i);
    A(i, :) = kernel_int(t, xgrid);
end
A = A*dx; 

x_true = sin(xgrid)+xgrid.^1;  
b_true = A*x_true; 
b     = b_true + randn(m,1)*norm(b_true)*dx*0.01; 


%% iDARR: iterative Data-Adaptive RKHS regularization
plotON   = 1; 
[x_reg,res,eta,k_corner] = iDARR_Lcurve(A,b,expl_measure,plotON); 


x_lse = lsqminnorm(A,b); 
x_pinv = pinv(A,1e-6)*b; % also very bad
figure; 
plot(xgrid, x_true,'LineWidth',1); hold on
plot(xgrid,x_reg,'--','LineWidth',1.5); 
plot(xgrid,x_lse,':','LineWidth',1.5); 
% plot(xgrid,x_pinv,':','LineWidth',1); 
ylim([min(x_true),max(x_true)]);
legend('True','iDARR','lsqmininorm'); 



%%-----------------------
function [x_reg,res,eta,k_corner] = iDARR_Lcurve(A,b,expl_measure,plotON)
% Exploration measure
if strcmp(expl_measure, 'L1')
   rho = sum(abs(A),1); 
elseif strcmp(expl_measure, 'L2')
   rho = sum(A.^2,1);
else
   error('Wrong exploration measure')
end

rho = rho'/sum(rho); 
figure; plot(rho); 

% IR with gGKB  
max_iter = 0.5*length(A(1,:)); 
reorth   = 1;
[X, res, eta] = dartr_spr(A, b, rho, max_iter, reorth);     % gGKB and updating procedure
[k_corner,info]   = Lcurve(res,eta,plotON,'l2-xHG');        % Early stopping by Lcurve
if isempty(k_corner)
   k_corner     =  find(log10(res)>-10,1,'last')-1; 
end
x_reg       = X(:,k_corner);

end