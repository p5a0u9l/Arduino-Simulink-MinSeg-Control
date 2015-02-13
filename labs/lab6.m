clear
%% Demo Problem 1
% a)
num = [0, 0, 1, 9, 26, 24];
den = [1, 14, 73, 176, 196, 80];

[A, B, C, D] = tf2ss(num, den);

% b) Asymptotically stable if Re{lambda_i} > 0 
if all(eig(A)) >= 0
    disp('System is Asymptotically stable')
else
    disp('System is Not Asymptotically stable')
end

% c), d)
CM = [];
for i = (1:size(A, 1)) - 1
    CM = [CM, A^i*B];
end

if rank(CM) >= rank(A)
    disp('System is CC')
else
    disp('System is not CC')
end

OM = [];
for i = (1:size(A, 1)) - 1
    OM = [OM; C*A^i];
end

if rank(OM) >= rank(A)
    disp('System is CO')
else
    disp('System is not CO')
end

% e)
h = tf(num, den);
sysr = minreal(h);
[num, den] = tfdata(sysr);
[Amin, Bmin, Cmin, Dmin] = tf2ss(num{:}, den{:});

% f)
Wc = gram(ss(sysr), 'c')
Wo = gram(ss(sysr), 'o')

% g)
x0 = [-50; 40; -300];
x1 = [0; 0; 0];
t1 = 5;
syms t
u = vpa(-Bmin.'*expm(Amin.'*(t1-t))*inv(Wc)*(expm(Amin*t1)*x0 - x1), 4);
% t = 0:0.01:t1;
u = eval(subs(u, t,0:0.01:t1)); 
% u = zeros(length(t), 1);
[y, t, x] = lsim(ss(sysr), u, 0:0.01:t1, x0);

subplot(121)
plot(t, u)
subplot(122)
plot(t, [y, x])
distFig

%% Individual Problem 
% a)
num = [0, 0, 0, 1, 3, 2];
den = [1, 15, 85, 225, 274, 120];

[A, B, C, D] = tf2ss(num, den);

% b) Asymptotically stable if Re{lambda_i} > 0 
if all(eig(A)) >= 0
    disp('System is Asymptotically stable')
else
    disp('System is Not Asymptotically stable')
end

% c), d)
cm = ctrb(A, B);
om = obsv(A, C);

if rank(cm) >= rank(A)
    disp('System is controllable')
else
    disp('System is not controllable')
end

if rank(om) >= rank(A)
    disp('System is observable')
else
    disp('System is not observable')
end

% e)
h = tf(num, den);
sysr = minreal(h);
[num, den] = tfdata(sysr);
[Amin, Bmin, Cmin, Dmin] = tf2ss(num{:}, den{:});

% f)
Wc = gram(ss(sysr), 'c')
Wo = gram(ss(sysr), 'o')

% g)
x0 = [10; 20; -30];
x1 = [0; 0; 0];
t1 = 5;
syms t
u = vpa(-Bmin.'*expm(Amin.'*(t1-t))*inv(Wc)*(expm(Amin*t1)*x0 - x1), 4);
% t = 0:0.01:t1;
u = eval(subs(u, t,0:0.01:t1)); 
% u = zeros(length(t), 1);
[y, t, x] = lsim(ss(sysr), u, 0:0.01:t1, x0);
figure, 
subplot(121)
plot(t, u)
subplot(122)
plot(t, [y, x])
distFig
