%% EE547 - Hw 4
% prepared by Paul Adams

function hw4()

%% Problem 1

    syms s t
    A = [14, -75, 190, -224, 96;
        1, 0, 0, 0, 0;
        0, 1, 0, 0, 0;
        0, 0, 1, 0, 0;
        0, 0, 0, 1, 0];
    
    % inverse Laplace Transform solution
    t1 = timeit(@() inverse_laplace_soln(A));
    Phi1 = inverse_laplace_soln(A)
    
    % matrix exponential solution
    t2 = timeit(@() expm(A*t));
    Phi2 = expm(A*t)
    
    % Jordan form solution
    [V, J] = jordan(A);
    Phi3 = V*expm(J*t)*inv(V)
    t3 = timeit(@() expm(J*t));
    
    % Function of a square matrix
    t4 = timeit(@() of_a_square_matrix(A));
    Phi4 = of_a_square_matrix(A)
    
    %% Timing results
    fprintf('Median Jordan form solution:                  %5.3f ms\n', t3*1000)
    fprintf('Median Exponential Matrix solution:           %5.3f ms\n', t2*1000)
    fprintf('Median Inverse Laplace solution:              %5.3f ms\n', t1*1000)
    fprintf('Median Functions of a square matrix solution: %5.3f ms\n', t4*1000)

%% Problem 2
    A = [-12 -55 -120 -124 -48;
        1 0 0 0 0;
        0 1 0 0 0;
        0 0 1 0 0;
        0 0 0 1 0];
    % a
    syms s t
    disp('Eigenvalues of A:')
    disp(eig(A))
    disp('Characteristic polynomial of A:')
    disp(det(s*eye(size(A, 1)) - A));
    % b
    [Q, J] = jordan(A);
    disp('Similarity matrix of A:')
%     Q = round(1e4*Q)/1e4;
    disp(Q)
    disp('Jordan form of A:')
    disp(J)
    At = A*t;
    Jt = J*t;
    f = (Q*expm(Jt))\Q + (Q*Jt)\Q;
    disp('f(A):')
    disp(f)
    

%% Functions
%% OF_A_SQUARE_MATRIX
function x = of_a_square_matrix(A)
eigvals = round(1e4*roots(charpoly(A)))/1e4;  
syms lambda t
%%
% compute $f$ on spectrum of $\mathbf{A}$
f = symfun(exp(lambda*t), lambda);
f_ = g_of_lambda(eigvals, f);
%%
% construct $h(\lambda) = \beta_0 + \beta_1\lambda + \beta_2\lambda^2 + \beta_3\lambda^3 + \beta_4\lambda^4$
h = symfun(lambda.^(0:4), lambda);
h_ = g_of_lambda(eigvals, h);
%%
% solve the system with $A*\beta = b$
% where $\beta$ are the coefficients of $h(\lambda)$
beta = h_\f_; 
%%
% solve for $h(\mathbf{A}) = f(\mathbf{A})$
x = zeros(5);
for i = 1:size(A, 1)
    x = x + beta(i)*A^(i-1);
end

%% G_OF_LAMBDA
function y = g_of_lambda(eigvals, g)
syms lambda
eigvals = sort(eigvals);
% get the multiplicies of eigvals
n = hist(eigvals, max(eigvals));
idx = 1; 
for i = 1:length(n)
    l = eigvals(i);
    % take derivative to the n-1 for eigenvalue with multiplicity n
    for j = 1:n(i)
        y(idx, :) = diff(g(lambda), lambda, j-1);
        y(idx, :) = subs(y(idx, :), lambda, eigvals(i));
        idx = idx + 1;
    end
end

%% INVERSE_LAPLACE_SOLN
function y = inverse_laplace_soln(A)
syms s t
sI_A = s*eye(size(A,1)) - A;
y = ilaplace(inv(sI_A));
