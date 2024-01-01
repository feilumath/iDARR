
% Using 2D large-scale test problems, and show the effect of the algorithms
% for reconstructing solutions.

clear, clc; close all; 
add_mypaths_discrete; 
rng(2023); 

% create image data
img = imread('toyobjects.png');  % toyobjects.png   'HSTgray.jpg'
img1 = imresize(img,[128,128]);
optblur.trueImage = img1;
% optblur.BlurLevel = 'mild';
optblur.CommitCrime = 'on';
optblur.BC = 'zero';
[A, b_true, x_true, ProbInfo] = PRblurspeckle(optblur);
% add noise
nel = 1e-2;  % noise level
b = PRnoise(b_true, 'gauss', nel);  % Observed Image

% prepare algorithms
[m, n] = sizm(A);
xn = norm(x_true);

rho = zeros(n,1);
for i = 1:n
    ei = zeros(n,1);
    ei(i) = 1;
    ai = mvp(A, ei);
    rho(i) = sum(abs(ai));
end
% rho = rho / sum(rho);

% compare reguarization methods
k = 200;
er1 = zeros(k,1);
er2 = zeros(k,1);
er3 = zeros(k,1);

% [X1, stop1] = idarr_deblur(A, b, k);
[X1, res1, eta1] = dartr_spr1(A, b, rho, k, 1);
[X2, res2, eta2] = lsqr_b(A, b, k, 1);
options3 = IRset('NoStop', 'on', 'RegParam', 'wgcv',...
    'NoiseLevel', nel, 'Reorth','on');
[X3, info3] = IRhybrid_lsqr(A, b, 1:k, options3);

[stop1, info1] = Lcurve(res1,eta1,1,'l2-xHG');
[stop2, info2] = Lcurve(res2,eta2,1,'l2-l2');
stop3 = info3.StopReg.It;

for i =1:k
    er1(i) = norm(x_true-X1(:,i)) / xn;
    er2(i) = norm(x_true-X2(:,i)) / xn;
    er3(i) = norm(x_true-X3(:,i)) / xn;
end

[~, k01] = min(er1);
[~, k02] = min(er2);
[~, k03] = min(er3);


%----------- plot -----------------------------------
figure;
PRshowx(x_true, ProbInfo)
title('True image','interpreter','latex','fontsize',18)
set(gca,'fontsize',18)

figure;
PRshowx(b, ProbInfo)
title('Noisy data','interpreter','latex','fontsize',18)
set(gca,'fontsize',18)

figure;
PRshowx(X1(:,k01), ProbInfo)
title(['Best solution, iDARR, $k_{0}$ = ',num2str(k01)],...
    'interpreter','latex','fontsize',18);
set(gca,'fontsize',18)

figure;
PRshowx(X2(:,k02), ProbInfo)
title(['Best solution, LSQR, $k_{0}$ = ',num2str(k02)],...
    'interpreter','latex','fontsize',18);
set(gca,'fontsize',18)

figure;
PRshowx(X1(:,stop1), ProbInfo)
title(['iDARR recons. image, $k$ = ',num2str(stop1)],...
    'interpreter','latex','fontsize',18);
set(gca,'fontsize',18)

figure;
PRshowx(X2(:,stop2), ProbInfo)
title(['LSQR recons. image, $k$ = ',num2str(stop2)],...
    'interpreter','latex','fontsize',18);
set(gca,'fontsize',18)

figure;
semilogy(1:k, er1, '->','Color','b','MarkerIndices',1:9:k,...
    'MarkerSize',6,'MarkerFaceColor','b','LineWidth',1.5);
hold on;
semilogy(1:k, er2, '-s','Color',[1,0.47,0.1],'MarkerIndices',1:9:k,...
    'MarkerSize',6,'MarkerFaceColor',[1.0,0.47,0.1],'LineWidth',1.5);
hold on;
semilogy(1:k, er3, '-o','Color','g','MarkerIndices',1:9:k,...
    'MarkerSize',6,'MarkerFaceColor','g','LineWidth',1.5);
plot(stop1, er1(stop1),'ko', 'MarkerSize',16, 'LineWidth',2)
plot(stop2, er2(stop2),'ko', 'MarkerSize',16, 'LineWidth',2)
plot(stop3, er3(stop3),'ko', 'MarkerSize',16, 'LineWidth',2)
xlabel('Iteration','fontsize',16);
legend('iDARR', 'LSQR','hybrid-l2', 'Location', 'northeast','fontsize',16);
ylabel('l2 relative  error','fontsize',16);
grid on;
grid minor;
set(gca, 'GridAlpha', 0.3);
set(gca, 'MinorGridAlpha', 0.01);
