clear
close all
%% demo
% 
% %%
% A = [1, 0, 0, 0; ...
%     2, 2, 0, 0; ...
%     0, 0, 3, 3; ...
%     0, 0, 0, 4;];
% 
% Adet = det(A);
% syms s t
% sI_A = s*eye(size(4)) - A;
% [V, D] = eig(A);
% x1 = expm(A*t);
% x2 = ilaplace(inv(sI_A));
% d = diag(expm(A*t));
% for i =1:length(d)
%     ezplot(d(i), [0 1]), hold on
% end
% hold off
% % distFig

%% individual
A = [0, -2; 2, -4];
syms t s
stm = expm(A*t);
CharMat = inv(s*eye(size(A)) - A);
CharPoly = poly(A);
[V, D] = eig(A);
d = diag(stm);
figure, hold on
for i = 1:length(d)
    ezplot(d(i), [0 5])
end
hold off

B = [1; -2];
C = eye(size(A));
sys = ss(A, B, C, 0);
impulse(sys)