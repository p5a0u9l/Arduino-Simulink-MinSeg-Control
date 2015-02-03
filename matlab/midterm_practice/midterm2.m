%% Problem 2
A1 = [1, 0; -3, 2];
B1 = [1; 0];
C1 = [1, 1];
D1 = 0;
A2 = [-2, -5; 0, -1];
B2 = [1; 0];
C2 = [1, 1];
D2 = 0;

% a) Zero-state equivalent?
% do not have same transfer function or eigenvalues, so... no.

%% Problem 3
A = [3, 5; 5, 3];
B = [1; 1];
C = [0, 1];
D = 0;
syms s t tau
Gs = C*inv(s*eye(2) - A)*B + D;
pretty(Gs)
t0 = 0;
Phi = expm(A*(t - t0));
disp(simplify(Phi))

zir = C*Phi*[1; 1]

zsr = C*int(expm(A*(t - tau))*B, tau, 0, t) 

pretty(zir + zsr)