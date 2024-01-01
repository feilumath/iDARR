function sysInfo = system_settings(kernel_type)
% basic settings of the weighted transformation

lb = 1;
rb = 5;      % Space Integration Interval
xn = 100;  


T = 5;      % Time integration interval
tn = 500;   % Time observation steps



% f_true_func = @(x) (sin(x-6)).^2-3;   % True f
% f_true_func = @(x) (x-6).^2+3;   % True f

% f_true_func = @(x) 15*((sin(x-6)).^2-3);   % True f
% f_true_func = @(x) 15*((sin(2*x-6)).^2-3);
% sigma = 0.001;    % noise level

%%
sysInfo.lb = lb;
sysInfo.rb = rb;
sysInfo.T = T;
sysInfo.tn = tn;
sysInfo.xn = xn;

switch kernel_type
    case 'exp'
        sysInfo.phi = @(t, x) x.^(-2).*exp(-x.*t); % This phi is the integral kernel K(t,x)
        sysInfo.kernel_type = 'exp';
    case 'poly'
        sysInfo.phi = @(t, x) x.^(-1).*abs(sin(x.*t +1)); % This phi is the integral kernel K(t,x)
        sysInfo.kernel_type = 'poly';
end


sysInfo = update_system_settings(sysInfo);


end