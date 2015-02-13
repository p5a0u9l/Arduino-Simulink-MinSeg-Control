%% Homework 5 - EE 547 (PMP) - Winter 2015
% prepared by Paul Adams
%

%% Initialization
function hw5()
close all
digits(3);
opengl('save', 'software')
format shortG
set(0, 'defaultTextInterpreter', 'latex'); 
numerical_precision = 1e-9;
syms s xdot x x_1 x_2 x_3 x_4 x_5 u_1 u_2 u_3 u_4 u_5 u t t_0


%% Problem 1
%%
% <html> <h3> Setup </h3> </html>
x = [x_1; x_2; x_3; x_4; x_5];
u = [u_1; u_2; u_3; u_4; u_5];
A = [-9 -31 -51 -40 -12;
     1 0 0 0 0;
     0 1 0 0 0;
     0 0 1 0 0;
     0 0 0 1 0];
B = [1; 0; 0; 0; 0];
C = [95 92 -3 60 -72];
D = 0;
render_latex(['\dot{x} = ' latex(sym(A)) '\mathbf{x} + ' latex(sym(B)) 'u'], 12, 1.3)
render_latex(['y = ' latex(sym(B)) '\mathbf{x} + ' latex(sym(D)) 'u'], 12, 1.3)
%%
% <html> <h3> a) Check if the system is asymptotically stable. </h3> </html>
%%
% System is asymptotically stable if $\Re{\{\lambda_i\}} > 0$ for all
% eigenvalues
if all(eig(A)) >= 0
    disp('System is Asymptotically stable')
else
    disp('System is Not Asymptotically stable')
end
%%
% <html> <h3> b) Evaluate the controllability and observability matrices of the system, C and O. </h3> </html>
cm = ctrb(A, B);
om = obsv(A, C);
render_latex(['\mathcal{C} = ' latex(sym(cm))], 12, 1.3)
render_latex(['\mathcal{O} = ' latex(sym(om))], 12, 1.3)
%%
% <html> <h3> c) Check if the system is controllable and observable by the ranks of controllability and observability matrices, C and O. </h3> </html>
n = size(A, 1);
if rank(cm) >= n
    disp('System is controllable')
else
    disp('System is not controllable')
end

if rank(om) >= n
    disp('System is observable')
else
    disp('System is not observable')
end
%%
% <html> <h3> d) Find the controllability and observability gramians, Wc and Wo, of this system.</h3> </html>
sys = ss(A, B, C, D);
Wc = gram(sys, 'c')
Wo = gram(sys, 'o')
%%
% <html> <h3> Please derive input that drives initial state x0 into x1 within 8 seconds.</h3> </html>
x0 = [-50; 40; -300; -100; 200];
x1 = zeros(5, 1);
t1 = 8;
syms t
u = vpa(-B'*expm(A'*(t1-t))*inv(Wc)*(expm(A*t1)*x0 - x1), 4);
render_latex(['u = ' latex(u)], 11, 0.5)
%%
% <html> <h3> Use lsim function to simulate this controllable system and plot the input and state variable. </h3> </html>
tspan = 0:0.01:t1;
u = eval(subs(u, t, tspan)); 
[y, t, x] = lsim(ss(sys), u, tspan, x0);

figure, 
plot(t, x)
grid on
title('State variables')
xlabel('time')
legend('x_1', 'x_2', 'x_3', 'x_4', 'x_5')

figure
grid on
plot(t, u)
title('Applied Input')
xlabel('time')

%% Problem 2
%%
% <html> <h3> Find Transfer Function for System 1 </h3> </html>
  