%% Lab2 Adams
%% iniz. 
close all
clear
titles = {'Step Input', 'Sine Input', 'Ramp Input'};
x_ini = [0.05; 0.05];
dt = 0.01;
t = 0:dt:10-dt;

%% A B C D
A = [0, 1; 0, -3];
B = [0; 1];
C = eye(2);
D = zeros(2, 1);

%%
[num, den] = ss2tf(A, B, C, D, 1);
[AA, BB, CC, DD] = tf2ss(num, den);

%%
tf_ = tf({num(1, :), num(2, :)}, den);
z = zero(tf_);
p = pole(tf_);

%%
input = 1;
sim('lab2S.slx', t(end));
%% plots
figure
plot(clk, [u, y])
title(titles{input})
legend({'Input', 'Output 1', 'Output 2'})

input = 2;
sim('lab2S.slx', t(end));
%% plots
figure
plot(clk, [u, y])
title(titles{input})
legend({'Input', 'Output 1', 'Output 2'})

input = 3;
sim('lab2S.slx', t(end));
%% plots
figure
plot(clk, [u, y])
title(titles{input})
legend({'Input', 'Output 1', 'Output 2'})

distFig('rows', 2, 'columns', 2)
