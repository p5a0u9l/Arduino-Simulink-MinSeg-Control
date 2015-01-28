clear 
close all

%% Problem 4
syms t
A = [-7, -12; 1, 0];
B = [1; 0];
C = [-4, -10];
D = 0;
x0 = [1; 2];

%% Solution to part a
pretty(expm(A*t)*x0)

%% Solution to part b
sys = ss(A, B, C, D);
impulse(sys)
