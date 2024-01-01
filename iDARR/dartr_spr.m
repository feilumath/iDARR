function [X, res, eta] = dartr_spr(A, b, rho, k, reorth,printON)
    % Subspace projection regularization method for adaptive RKHS norm regularization.
    % The regularizer is constructed by kernel dependent RKHS as
    %       ||x||_{G}^2 with G^{-1} = B^{-1}A^TAB^{-1},
    % where B is a diagonal matrix whose elements are the atom of
    % exploration measure, which can be computed as
    %      B(i,i) = rho(i).
    % It computes iterative solutions of
    % min||x||_G  s.t. min||Ax-b||_2,   G^{-1} = B^{-1}A^TAB^{-1}.
    % The Cholesky G=L^TL need not be computed explicity.
    % It is mathematically equivalent to LSQR applied to
    %   ||b - AL^{-1}(Lx)||_2
    % by the Lanczos bidiagonaliation of AL^{-1}.
    %
    % Inputs:
    %   A: either (a) a full or sparse mxn matrix;
    %             (b) a matrix object that performs the matrix*vector operation
    %  rho: atoms of the exploration measure
    %   b: right-hand side vector
    %   A and b construct the ill-posed linear system: Ax + e = b, where e is the noise
    %   k: the maximum number of iterations  
    %   reorth: 
    %       0: no reorthogonalization
    %       1: full reorthogonaliation, MGS
    %       2: double reorthogonaliation, MGS
    %
    % Outputs: 
    %   X: store the first k regularized solution
    %   res: strore residual norm of the first k regularized solution
    %   eta: strore solution G-norm of the first k regularized solution
    %
    % Reference: [1]. Fei Lu, M-J Y. Ou, An adaptive RKHS regularization for 
    %   Fredholm integral equations, arXiv:2303.13737v2, 2023.
    % [2]. Haibo Li, A preconditioned Krylov subspace method for linear inverse 
    %   problems with general-form Tikhonov regularization, preprint, 2023.
    %
    % Haibo Li, 04, Aug, 2023.
    % 


    % Check for acceptable number of input arguments
    if nargin < 5
        error('Not Enough Inputs')
    end
    if ~exist("printON",'var'); printON = 0; end

    rho = rho(:);
    % A2 = A' * A;  % if the matrix form of A is not available, A2 can not be explicitly constructed
    [m, n] = size(A); 
    if m ~= size(b,1) || n ~= size(rho,1)
        error('The dimensions are not consistent')
    end
      
    tol_alpha_beta = 1e-9; %% set the tolerance of small alpha, beta. Avoid using small eigenvalues's eig-spaces.  
    k_lowbound = 5;         % the lower bound of k dimensions to be computed, so that the L-curve will work.  
    tol_res    = 1e-15; 
    tol_eta    = 1e-6; 

    % declares the matrix size
    B = zeros(k+1, k);
    U = zeros(m, k+1);
    Z = zeros(n, k+1);
    P = zeros(n, k+1);  % P = GZ, assisting with reorthogonalization of z (avoid G*z computations)

    X = zeros(n, k);
    res = zeros(k ,1);  % residual norm
    eta = zeros(k ,1);  % solution norm
    
    % start iteration, compute u_1, z_1
    % fprintf('Start the Subspace projection iteration ============\n');
    bbeta = norm(b);
    u = b / bbeta;  U(:,1) = u;
    r = A' * u;
    r1 = A' * (A * (r./rho));
    rz = r1 ./ rho;
    % rz = (A2 * (r./rho)) ./ rho;
    alpha = sqrt(r' * rz);  B(1,1) = alpha;
    z = rz / alpha;  Z(:,1) = z;
    p = r / alpha;   P(:,1) = p;

    % Prepare for update procedure
    w = z;
    phi_bar = bbeta;
    rho_bar = alpha;
    x = zeros(n, 1);
    x_bar  = zeros(n, 1);  % x_bar = Gx, assist to compute ||x||_G
    w_bar = p;  % w_bar = Gw, used to iteratively update x_bar

    beta = 1; j=0; 
    % The j-th step gen-GKB iteration and update procedure
    while  j<=k && alpha*beta >0  
        j= j+1; 
        % compute u in 2-inner product
        if printON==1
            fprintf('Running the SPR process: the %d-th step ---------\n', j);
        end
        s = A * z - alpha * u;
        if reorth == 1  % full reorthogonalization of u, in 2-inner product
            for i = 1:j 
                s = s - U(:,i)*(U(:,i)'*s);
            end
        elseif reorth == 2  % double reorthogonalization of u
            for i = 1:j 
                s = s - U(:,i)*(U(:,i)'*s);
            end
            for i = 1:j 
                s = s - U(:,i)*(U(:,i)'*s);
            end
        else
            % pass
        end
        beta = norm(s);  B(j+1,j) = beta;                    
        if j>k_lowbound && (isnan(beta) || beta <= tol_alpha_beta); j= j-1;  break; end
        u = s / beta;    U(:,j+1) = u;

        r = A' * u - beta * p;
        if reorth == 1  % full reorthogonalization of p over z_i
            for i = 1:j 
                r = r - Z(:,i)*(Z(:,i)'*r);
            end
        elseif reorth == 2
            for i = 1:j 
                r = r - Z(:,i)*(Z(:,i)'*r);
            end
            for i = 1:j 
                r = r - Z(:,i)*(Z(:,i)'*r);
            end
        else
            % pass
        end
        r1 = A' * (A * (r./rho));
        rz = r1 ./ rho;
        % rz = (A2 * (r./rho)) ./ rho;
        alpha = sqrt(r' * rz);   B(j+1,j+1) = alpha;          
        if j>k_lowbound && (isnan(alpha)|| alpha<= tol_alpha_beta); j= j-1; break;  end
        z = rz / alpha;  Z(:,j+1) = z;   
        p = r / alpha;   P(:,j+1) = p;

        % Construct and apply orthogonal transformation
        rrho = sqrt(rho_bar^2 + beta^2); 
        c1 = rho_bar / rrho;
        s1 = beta / rrho; 
        theta = s1 * alpha; 
        rho_bar = -c1 * alpha;
        phi = c1 * phi_bar;
        phi_bar = s1 * phi_bar;  

        % Update the solution and w_i
        x = x + (phi/rrho) * w;  
        w = z - (theta/rrho) * w;
        X(:,j) = x;
        res(j) = abs(phi_bar);  % ||b-A*x||_2

        x_bar = x_bar + (phi/rrho) * w_bar;  % G*x
        w_bar = p - (theta/rrho) * w_bar;    % G*w
        eta(j) = sqrt(x'*x_bar);
        if (j>k_lowbound) && (res(j)< tol_res || eta(j)-eta(j-1)<tol_eta); break;  end
    end
    
    
n_terminate = j; 
X   = X(:,1:n_terminate); 
res = res(1:n_terminate); 
eta = eta(1:n_terminate);
end
    
    
