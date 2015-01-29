clear 
close all

%% Problem 4
syms t
A = [-7, -12; 1, 0];
B = [1; 0];
C = [-4, -10];
D = 1;
x0 = [1; 2];

%% Solution to part a
Phi = expm(A*t);

%% Solution to part b
sys = ss(A, B, C, D);
impulse(sys)
