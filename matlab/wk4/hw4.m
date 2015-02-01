function hw4()
%% Problem 1
%%
% 
% <<p1.PNG>>
% 

syms t
A = [14, -75, 190, -224, 96;
    1, 0, 0, 0, 0;
    0, 1, 0, 0, 0;
    0, 0, 1, 0, 0;
    0, 0, 0, 1, 0];

% inverse Laplace Transform solution
tic
Phi = inverse_laplace_soln(A); toc
disp(Phi)

% matrix exponential solution
tic
Phi = expm(A*t); toc
disp(simplify(Phi))

%% Functions
function y = inverse_laplace_soln(A)
syms s
sI_A = s*eye(size(A,1)) - A;
y = ilaplace(sI_A);
