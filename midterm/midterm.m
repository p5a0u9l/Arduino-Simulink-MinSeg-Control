%% EE547 (PMP) Midterm - Winter 2015
% prepared by Paul Adams
%

%% Initialization
function midterm()
close all
opengl('save', 'software')
format shortG
set(0, 'defaultTextInterpreter', 'latex'); 
numerical_precision = 1e-9;
syms s x x_1 x_2 x_3 u t t_0

%% Problem 1
%%
% <html> <h3> a) Linearize the system around equilibrium points </h3> </html>
%
x = [x_1; x_2; x_3];
f = [-9*x_1 - 4*x_2 - (1+x_3)*x_3 + sin(x_3) + sin(u);...
    (x_2*x_3 - 4)*x_1 - 10*sin(x_2) + 3*cos(x_3) + x_3^2*sin(u);...
    9*x_1 + (x_1^2 - 4)*x_3 - 10*x_2 + u];
g = [x_1 + x_2*x_3 + sin(u);...
    x_2 + x_1*x_3 + u^2;...
    x_3 + x_2*x_3 + cos(u)];
render_latex(['f = ' latex(f)], 12, 0.7)
render_latex(['g = ' latex(g)], 12, 0.7)
%%
% The state-space matrices are found using 
% 
% $$\mathbf{A} = \left. \frac{\partial{f}}{\partial{\mathbf{x}}}\right|_{x^{eq}, u^{eq}}
% \mathbf{B} = \left. \frac{\partial{f}}{\partial{u}}\right|_{x^{eq}, u^{eq}}
% \mathbf{C} = \left. \frac{\partial{g}}{\partial{\mathbf{x}}}\right|_{x^{eq}, u^{eq}}
% \mathbf{D} = \left. \frac{\partial{g}}{\partial{u}}\right|_{x^{eq}, u^{eq}}
% $$
% 
% where, 
%%
% $$\mathbf{x} = \left[\begin{array}{c} x_{1}\\ x_{2}\\ x_{3} \end{array}\right]$$
%
%%
xeq = [-0.1 0.1 -0.2];
ueq = 0;
A = subs(jacobian(f, x), [x_1, x_2, x_3, u], [xeq, ueq]);
B = subs(jacobian(f, u), [x_1, x_2, x_3, u], [xeq, ueq]);
C = subs(jacobian(g, x), [x_1, x_2, x_3, u], [xeq, ueq]);
D = subs(jacobian(g, u), [x_1, x_2, x_3, u], [xeq, ueq]);
render_latex(['\mathbf{A} = ' latex(simplify(A))], 12, 0.75)
render_latex(['\mathbf{B} = ' latex(simplify(B))], 12, 0.75)
render_latex(['\mathbf{C} = ' latex(simplify(C))], 12, 0.75)
render_latex(['\mathbf{D} = ' latex(simplify(D))], 12, 0.75)
A = double(A);
B = double(B);
C = double(C);
D = double(D);
%%
% <html> <h3> b) Find the transfer functions of this system. 
% Are these transfer functions proper rational functions? </h3> </html>
%
% Given the state-space matrices $\mathbf{A}, \mathbf{B}, \mathbf{C}$ and 
% $\mathbf{D}$. The transer function is found using 
%
% $$\hat{G}(s) = \mathbf{C}(s\mathbf{I} - \mathbf{A})^{-1}\mathbf{B} $$
% 
% *NOTE* since the symbolic representation of G using 
% 
%   sI_A = s*eye(size(A)) - A;
%   G = C*inv(sI_A)*B;
% 
% results in unreadable expressions, use _ss2tf_ instead. 
[num, den] = ss2tf(A, B, C, D)
%%
% *Proper rational?* 
%   The degree of the numerators are, at most, equal to the degree of the
%   denominator. Therefore, the transfer functions are proper rational.
%%
% <html> <h3> c) Is the system BIBO stable?</h3> </html>
%
% As noted above, the system is proper rational. 
poles = roots(charpoly(A));     % eig is preferable, but is used below... 
if all(real(poles) < 0)
    fprintf('Since the poles (eigenvalues) each have negative real parts, \n the system is BIBO stable.\n')
else
    disp('The system is NOT BIBO stable.')
end
%%
% <html> <h3> d) Evaluate the Jordan form of the matrix A 
% and corresponding transformation matrix Q.</h3> </html>
%
%% 
% Find the eignevalues of A
%%
% *NOTE* use _eig_ to get more accurate results
lambda = eig(A)
%% 
% Since the eigenvalues of A are distinct, the Jordan form is simply the
% eigenvalues of A along the diagonal. 
% The tranformation matrix, Q, is found by finding a solution to the
% homogenuous equation
%%
% 
% $\left (\mathbf{A} - \lambda_i \mathbf{I} \right ) q_i = \mathbf{0}$
% 
J = diag(lambda);
Q = zeros(size(A));
for i = 1:length(lambda)
    Q(:, i) = null(A - lambda(i)*eye(size(A)));
end
J_ = array2table(J)
Q_ = array2table(Q)
%%
% <html> <h3> e) Find the state transition matrix from the Jordan form. </h3> </html>
%
% For an LTI system, the State Transition Matrix $\Phi(t, t_0)$ using the
% Jordan form is given by
%%
% 
% $$\Phi(t, t_0) = \mathbf{Q}e^{\mathbf{J}(t-t_0)}\mathbf{Q}^{-1}$$
% 
Phi = Q*expm(J*(t))*inv(Q);
Phi = vpa(Phi, 2)
% render_latex(['\Phi(t, 0) = ' latex(Phi)], 12, 0.75)
%%
% <html> <h3> f) Given the initial state evaluate and plot zero-input responses of
% the system. </h3> </html>
%
x0 = [-0.1; 0.1; -0.2];
sys = ss(A, B, C, D);
[y_zir, t, x_zir] = initial(sys, x0, 3);
%% 
% <html> <h3> g) Evaluate the complete impulse response of this system. Check 
% if the system is asymptotically stable by solving the Lyapunov equation, where
% Q is a symmetric positive-definite matrix: </h3> </html>
%
Q = eye(size(A, 1));
% The system is asymptotically stable if the solution to the Lyapunov
% equation, P, is positive definite. 
P = lyap(A, Q)
% Choose the condition that every eigenvalue of P be positive as the
% condition for the positive-definiteness of P
if all(eig(P) > 0)
    disp('P is positive definite, therefore, the system is asymptotically stable')
else
    disp('The system is NOT asymptotically stable')
end
% compute the impulse response
[y_zsr, ~, x_zsr] = impulse(sys, 3);
%%
% Plots
figure, grid on
initial(sys(1), sys(2), sys(3), x0, 3);
figure, grid on
impulse(sys(1), sys(2), sys(3), 3);
figure, grid on
plot(t, y_zir + y_zsr)
title('Complete Response')
xlabel('Time (seconds)'); ylabel('Amplitude')

%% Problem 2
%%
% <html> <h3> Find Transfer Function for System 1 </h3> </html>
%%
% Given the state-space matrices $\mathbf{A}, \mathbf{B}, \mathbf{C}$ and 
% $\mathbf{D}$. The transer function is found using 
%
% $$\hat{G}(s) = \mathbf{C}(s\mathbf{I} - \mathbf{A})^{-1}\mathbf{B} $$
% 
A = [-5 -9 4; 2 -9 2; 9 -10 -8];
B = [-1; 2; -3];
C = eye(size(A));
D = zeros(3, 1);
sI_A = s*eye(size(A)) - A;
G = C*inv(sI_A)*B;
render_latex(['\hat{G_1}(s) = ' latex(simplify(G))], 16, 1.2)
%%
% <html> <h3> Find Transfer Function for System 2 </h3> </html>
%
A_ = [-82 -8 54; -174 -33 130.5; -138 -16 93];
B_ = [2; -7; 0];
C_ = [-4 -1 3.5; 1 0 -0.5; 2 1 -2];
D_ = zeros(3, 1);
sI_A = s*eye(size(A_)) - A_;
G_ = C_*inv(sI_A)*B_;
render_latex(['\hat{G_2}(s) = ' latex(simplify(G_))], 16, 1.2)
%%
% <html> <h3> Verify equality using <i>ss2tf</i> </h3> </html>
%
G = ss2tf(A, B, C, D);
G_ = ss2tf(A_, B_, C_, D_);
if max(max(G - G_)) > numerical_precision
    disp('Systems are not zero-state equivalent')
else
    fprintf('Systems realize the same transfer function, therefore, \n they are zero-state equivalent.\n')
end

%% Problem 3
A = [-4 3 -1 10 3;...
    -6 -9 5 7 -6;...
    -8 -5 2 -9 2;...
    -4 -5 6 -1 2;...
    -6 8 -8 -4 2];
%%
% <html> <h3> a) Evaluate the eigenvalues and characteristic polynomial of <b>A</b> </h3> </html>
%
%%
% The characteristic polynomial of $\mathbf{A}$ is given by
%
% $$\Delta(\lambda) = \det{s\mathbf{I} - \mathbf{A}} $$
%
% And the eigenvalues of $\mathbf{A}$ are the roots of the characteristic
% polynomial
%
sI_A = s*eye(size(A)) - A;
CharPoly = det(sI_A);
render_latex(['\Delta(\lambda) = ' latex(CharPoly)], 12, 0.5)
lambda = roots(sym2poly(CharPoly))
%%
% <html> <h3> Verify by using <i>eig</i> and displaying </h3> </html>
%
scatter(real(lambda), imag(lambda)); hold on
scatter(real(eig(A)), imag(eig(A))); hold off
xlabel('Real Axis'); ylabel('Imaginary Axis'); 
title('Eigenvalues of matrix $\mathbf{A}$')
%%
% <html> <h3> b) Evaluate Jordon form and transformation matrix <b>Q</b> of matrix <b>A</b> </h3> </html>
%
% Since the eigenvalues of A are distinct, the Jordan form is simply the
% eigenvalues of A along the diagonal. 
% The tranformation matrix, Q, is found by finding a solution to the
% homogenuous equation
%%
% 
% $\left (\mathbf{A} - \lambda_i \mathbf{I} \right ) q_i = \mathbf{0}$
% 
lambda = eig(A);
J = diag(lambda);
Q = zeros(size(A));
for i = 1:length(lambda)
    Q(:, i) = null(A - lambda(i)*eye(size(A)));
end
J_ = array2table(J)
Q_ = array2table(Q)
%%
% <html> <h3> Verify transformation </h3> </html>
%
if max(max(A*Q - Q*J)) > numerical_precision
    disp('The transformation is not valid.')
else
    disp('The transformation is valid.')
end
%%
% <html> <h3> c) Evaluate function </h3> </html>
%
%%
% 
% $f_2(A) = A^5+10A^4+251A^3+1658A^2+12462^A+23160 I_{5x5}$ 
%
% The solution equates $f(\mathbf{A}) = f(\lambda)$ to $h(\lambda)$ where
% $h(\lambda) = \beta_0+\beta_1\lambda+...+\beta_{n-1}\lambda^{n-1}$
%
% However, I discovered the hard way that solving linear equations with 
% this $\mathbf{A}$ given the numeric precision available in Matlab leads to badly scaled
% matrices which factor with error. 
%
% Instead, use Jordan form of $\mathbf{A}$, then $\mathbf{A}^k = \mathbf{Q}^{-1} \mathbf{J}^k 
% \mathbf{Q}$. 
% 
% Since the eigenvalues of $\mathbf{A}$ are distinct, $\mathbf{J}^k = \lambda^k \mathbf{I_5}$, where $k = 5, 4, 3, 2, 1, 0$
% 
k = 5:-1:0;
b = [1 10 251 1658 12462 23160];
f2 = zeros(size(A));
f2_ = zeros(size(A));
for i = 1:length(k)
    A_to_the_k = inv(Q)*(diag(lambda.^k(i)))*Q;
    f2 = f2 + A_to_the_k*b(i);
    f2_ = f2_ + A^k(i)*b(i); % the direct solution
end
array2table(f2)
%%
% <html> <h3> Verify Results using explicit A^k </h3> </html>
%
if max(max(f2_ - f2)) > numerical_precision
    disp('The transformation is not valid.')
else
    disp('The transformation is valid.')
end