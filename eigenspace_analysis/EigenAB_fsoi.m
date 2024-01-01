function [V_A,eigA,V_L, eigL,r]= EigenAB_fsoi(A,B, plotON)
%% plot the eigenvalues of A and generalized (A,B)
method = 'Eig'; % 'SVD': 
    %  using eig(A) leads to similar V, but SVD ensures >= eigenvalues; 
    %  eig(A,B) is the right one to get V'*B*V = I 
switch method
    case 'Eig'
        [V_A, eigA]  = eig(A);      % A*V = V*D.     V'*V = I 
        [V_L, eigL0] = eig(A, B);   % A*V = B*V*D.   V'*B*V = I 
        [eigA, indA]  = sort(diag(eigA),'descend');
        [eigL, indAB] = sort(diag(eigL0),'descend');     
        V_A  = V_A(:, indA);
        V_L = V_L(:, indAB);
    case 'SVD'    % A = U*S*V';  U is "same" as V_A (different sign); but S is nonnegative
        [U,S,V]      = svd(A); eigA = diag(S); V_A = U;  
        [V_L, eigL0] = eig(A, B);   % A*V = B*V*D.   V'*B*V = I   
        [eigL, indAB] = sort(diag(eigL0),'descend');   V_L = V_L(:, indAB);
end

% r = rank(A);  % the threshold is 1e-16
r = length(find(abs(eigL)>1e-8));
if plotON ==1
    figure;
    subplot(131)
    semilogy(abs(eigA),'r-x','linewidth',1); hold on;
    semilogy(abs(eigL),'b-o','linewidth',1);
    legend('Eig val A','Eig val (A,B)');title('eigenvalues of A and (A,B)')
    
    
    subplot(132)
    plot(V_A(:,1:min(r,6)),'linewidth',1); title('eig vect A')
    subplot(133)
    plot(V_L(:, 1:min(r,6)),'linewidth',1); title('eig vect (A,B)')
end
end 


function compuareSVD_eig(A,B)
% compare SVD with eig, just checking the numerical difference 
[U,S,V]      = svd(A); %  diagonal matrix S, with nonnegative diagonal elements in decreasing order, and unitary matrices U and V so that A = U*S*V'.
S = diag(S); 

[V_A, eigA]  = eig(A); 
[eigA, indA] = sort(diag(eigA),'descend');
V_A          = V_A(:, indA); 

    figure;
    subplot(131)
    semilogy(abs(eigA),'r-x','linewidth',1); hold on;
    semilogy(abs(S),'b-o','linewidth',1);      legend('Eig','SVD'); 
    
    subplot(132)  % the eigen-vectors are the same, but with different signs 
    plot(V_A(:,1:3),'linewidth',1); title('eig vect A')
    subplot(133)
    plot(-U(:, 1:3),'linewidth',1); title('eig vect SVD')
    
    
%    subplot(132);  imagesc(V_A); title('V in Eig'); subplot(133); imagesc(U); title('U in SVD'); 
    
[V_L, eigL0]  = eig(A, B);
[eigL, indAB] = sort(diag(eigL0),'descend'); 
V_L           = V_L(:, indAB);


[U_L, S_L, V_svd]  = svd(B\A);   % here B is diagnal, so B\A is sysmetric 
figure; subplot(131); semilogy(abs(eigL),'r-x','linewidth',1); hold on;
        semilogy(S_L,'b-o','linewidth',1);      legend('EigAB','svd((invB)A)');
        subplot(132); plot(V_L(:,1:3),'linewidth',1); title('eig(A,B)')
        subplot(133); plot(U_L(:, 1:3),'linewidth',1); title('svd((invB)A)')

    
end


function compuareGSVD_Geig(A,B)
% we should use eig(A,B) to get V'*B*V = I
[U,V,X,C,S] = gsvd(A,B); % --- it is not eig(A,B) % unitary matrices U and V, a (usually) square matrix X, and nonnegative diagonal matrices
%   C and S so that 
%       A = U*C*X'
%       B = V*S*X'
%       C'*C + S'*S = I


[V_L, eigL0]  = eig(A, B); 
[eigL, indAB] = sort(diag(eigL0),'descend'); 
V_L           = V_L(:, indAB);


invBA = B\A; [V_AinvB, eigAinvB]  = eig(invBA);  [eigAinvB, indAB] = sort(diag(eigAinvB),'descend'); V_AinvB = V_AinvB(:, indAB);

[U, svdAB,V_svdAB] = svd(invBA); 
figure; subplot(141); semilogy(abs(eigL),'r-x','linewidth',1); hold on;
        semilogy(svdAB,'b-o','linewidth',1);     
        semilogy(eigAinvB,':d','linewidth',1.5);
        legend('EigAB','svd(invBA)','eig inBA');
        
        subplot(142); plot(V_L(:,1:3),'linewidth',1); title('eig(A,B)')
        subplot(143); plot(U(:, 1:3),'linewidth',1); title('svd(invBA)')
        subplot(144); plot(V_AinvB(:, 1:3),'linewidth',1); title('eig invBA)')
end

