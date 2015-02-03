%% Problem 2
A = [0, 3;3, 0];
B = [0; 1];
C = [1, 1];

% a) transfer function
sys = ss2tf(A, B, C, 0);
disp(sys)
% b) state transition matrix
t0 = 0;
syms t
Phi = expm(A*(t - t0));
disp(Phi)
% c) ZIR of the system
x0 = [1/2; 1/2];
x = Phi*x0;
disp(x)

%% Problem 3
A = [1, 0, 1; 0, 5, 2; 0, 0, 3];

% a) characteristic polynomial
charPoly = poly(A);
disp(charPoly)
% b) eigenvalues
lambda = round(roots(charPoly));
disp(lambda)
% c) Matrix power
k = 25;
% f(lambda) = lambda^k
f = lambda.^k;
% h(lambda) = b0 + b1*lambda + b2*lambda^2
for i = 1:length(lambda)
    h(i, :) = lambda(i).^(0:2);
end
beta = h\f;
Ak = beta(1)*eye(size(A, 1)) + beta(2)*A + beta(3)*A^2;
disp(Ak)
A^k - Ak;

%% Problem 4
A1 = [2, 6, 2; 2, 0, 0; 2, 6, 2];
A2 = [-1, 7, 4, 9; 0, 0, -2, -3; 0, 0, 0, 0; 0, 0, 0, -2];
% a) charpoly and eigenvalues
l1 = round(roots(poly(A1)));
l2 = round(roots(poly(A2)));
disp(l1)
disp(l2)
% b) Jordan form
% using rank(Ai - diag(li)), A1 is full rank, A2 has
disp(jordan(A1))
disp(jordan(A2))