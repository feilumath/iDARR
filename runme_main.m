  
%{
weighted Deconvolution in the form of inversion
                L f = y,   size(L) =  n_y x n_u       size(f) = n_u 
      << == >>  Ax  = b
Solution: least square with RKHS-regularization 
Key: space of identifiability, exploration measure 
%}


clc; close all; clear all;
add_mypaths_discrete;    % get SAVE_DIR = local dir for saving data
rng(1)


%% Load system settings
kernel_type = 'exp';    % choose the spectrum decay type: 'exp' or 'poly'
sysInfo    = system_settings(kernel_type);
A     = sysInfo.L_operator;
xn    = sysInfo.xn;
tn    = sysInfo.tn;
tgrid = sysInfo.tgrid;
xgrid = sysInfo.xgrid;
dx    = sysInfo.dx;

method = 'IR';    % choose the regularization method: 'IR' or 'Direct'

%% Get rho and L2(rho) basis matrix B 
% expl_measure = 'L1';  % exploration measure
expl_measure = 'L2';  % exploration measure
if strcmp(expl_measure, 'L1')
      rho = sum(abs(A),1); 
elseif strcmp(expl_measure, 'L2')
      rho = sum(A.^2,1);
else
      error('Wrong exploration measure')
end

rho = rho/(sum(rho)*dx);  % normalize, does not seem necessary for this example, maybe other examples
figure; % plot the exploration measure 
plot(xgrid, rho,'linewidth',1); xlabel('u');ylabel('rho');
B   = diag(rho);

%% Analysis function space of identifability
% It is for analysis purpose only. 
% Get normal matrix Abar; to get vector bbar later for different f
Abar                     = A'*A;  
[V_A,eigA,V_AB, eigAB,r] = EigenAB_fsoi(Abar,B,1); 

%% unconstained LSE with regularuzations: l2, L2, RKHS
% includes a single test for demonstration and tuning and multiple tests for robustness
reguLSE_unconstrained; 





