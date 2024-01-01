%% test the computational time for the IR-DA-RKHS, hybrid-RKHS, and DARTR.
% Roughly: O(k_iter *n^2), O(n^2*k^2), O(n^2)

clc; close all; clear all;
add_mypaths_discrete;    % get SAVE_DIR = local dir for saving data
rng(1)


%% test settings
xn_seq    =  [200,400,800,1600,3200]; %  [6400,12800];%
methods = {'Direct','IR'};

DnormType  = {'RKHS'};  
IRnormType = {'iDARR'};

kernel_type = 'poly'; 

nsimu     = 10; % number of simulations
n_methods = length(methods); 
n_xnseq   = length(xn_seq);

plotON = 0;
% filename =  [SAVE_DIR,'time_complexity_xnseq',num2str(n_xnseq),'.mat']; 

% if ~exist(filename,'file')
    sysInfo    = system_settings(kernel_type);
    tn         = sysInfo.tn;
    timetable  = zeros(nsimu,n_methods,length(xn_seq));
    for i = 1:n_xnseq
        sysInfo.xn = xn_seq(i);
        sysInfo = update_system_settings(sysInfo);
        A = sysInfo.L_operator;
        tgrid = sysInfo.tgrid;
        xgrid = sysInfo.xgrid;
        dx    = sysInfo.dx;
        %% Get regression matrix A,
        % to get vector b later for different f
        Abar = A'*A;

        %% Get rho and L2(rho) basis matrix B
        rho = sum(A);  rho = rho/(sum(rho)*dx);  % normalize, does not seem necessary for this example, maybe other examples
        B   = diag(rho);

        %% analysis function space of identifability
        [V_A,eigA,V_AB, eigAB,r]= EigenAB_fsoi(Abar,B,0);

        %% f_true outside the x) (sin(x-6)).^2+1;   % True f
        % f_true_func = @(x) 0.7*exp(-(x-2).^2/0.25)*sqrt(1/(2*.5*pi)) +.3*exp(-(x-4).^2/0.09)*sqrt(1/(2*0.3*pi));   % True f
        % f_true_func = @(x) 15*((sin(x-6)).^2-3);   % True f
        % f_true_func = @(x) 15*((sin(2*x-6)).^2-3);
        f_true_func = @(x) x.^2;
        f_true = f_true_func(xgrid);

        %%  each method simulate nsimu times
        time_temp = zeros(nsimu,n_methods);
        Lrep      = sysInfo.dt/sysInfo.T;
        parfor j = 1:nsimu
            nsr = 0.5;
            b_true = A*f_true;
            b_true_norm = sqrt(sum(b_true.^2)*Lrep);
            b      = b_true + b_true_norm*nsr*randn(tn, 1);

            for k = 1:n_methods
                method = methods{k};
                switch method
                    case 'Direct'
                        tic
                        normType = DnormType;
                        bbar      = A'*b;
                        [~,~,cLcurveAll] = reg_Lcurve_3in1(Abar,bbar,B,plotON,normType);

                    case 'IR'
                        tic
                        normType = IRnormType;  % iterative method
                        [~,cLcurveAll] = IRreg_3in1(A,b,B,rho,plotON,normType);
                end
                t1  =toc;
                time_temp(j,k)  = t1;
            end
            fprintf('(j-simu,xn)= (%d,%d) out of (%d,%d)  \n ', j, i,nsimu,n_xnseq);
        end
        timetable(:,:,i)  = time_temp;
        % save(filename,"timetable","xn_seq"); 
    end
% end

% load(filename); 

close all
newfigure = 1;
lgnd      = methods;
label_y   = 'Computational Time (seconds)'; 
fig_dir   = [SAVE_DIR,'figures/']; 
fig_name  = [fig_dir,'Computational_time']; 
data_array =  timetable;   % n_samples x n_types x n_ind
label_x    = 'n';
plot_scale = 'log';  % 'log' or 'linear'
plot_mean_std(xn_seq,data_array,newfigure,lgnd,label_y,label_x,plot_scale,fig_name);




