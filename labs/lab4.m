%% demo
echo on
A = [6 -12 8;
    1 0 0;
    0 1 0];
% compute the characteristic polynomial and eigenvalues of matrix A
syms lambda
n = rank(A);
char_poly = poly(A);
rts = roots(char_poly);
[x, lambda] = eig(A);
[Q, JA] = jordan(A);
Q\A*Q == Q\A*Q;

%% individual
x0 = [1 0 1 0 1]';
A = [14 -75 190 -224 96;
    1 0 0 0 0;
    0 1 0 0 0;
    0 0 1 0 0;
    0 0 0 1 0]; 
char_poly = poly(A);
[x, lambda] = eig(A);
[Q, JA] = jordan(A);

syms t
x = Q*expm(JA*t)\Q*x0;
hold on

for i=1:length(x)
    ezplot(x(i), [0 10])
    names{i} = char(x(i));
end
legend(names)
title('Five plots')
axis tight
hold off
x_check = expm(A*t)*x0;