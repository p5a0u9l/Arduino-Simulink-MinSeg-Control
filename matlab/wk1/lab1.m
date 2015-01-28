function lab1()

switch prob
    case 'demo'
        run_demo_prob
    case 'individual'
end

function run_demo_prob()
%% iniz.
t = 0:0.01:5;
xini = [0.5; 0.5];

%% A, B, C, D matrices
A = [0 1; 1 1];
B = [0; 0];
C = eye(2);
D = [0; 0];

%% create input
u = cos(t);

%% StateSpace Model
ss1 = ss(A, B, C, D);

%% Simulation
ZeroInputResponse = initial(ss1, xini, t);

[ImpulseInputResponse, t1] = impulse(ss1);
[StepInputResponse, te] = step(ss1);

LinearResponse = lsim(ss1, u, t, xini);

%% Plots
figure
subplot(221), plot(t, ZeroInputResponse)
xlabel('time [seconds]')
title('Zero Input Response')

subplot(222), plot(t1, ImpulseInputResponse)
xlabel('time [seconds]')
title('Impulse Input Response')

subplot(223), plot(te, StepInputResponse)
xlabel('time [seconds]')
title('Step Input Response')

subplot(224), plot(t, LinearResponse)
xlabel('time [seconds]')
title('Cosine Input Response')

function individual_problem
%% iniz.
t = 0:0.01:5;
xini = [0.5; 0.5];

%% A, B, C, D matrices
A = [0 1; 1 1];
B = [0; 0];
C = eye(2);
D = [0; 0];

%% create input
u = cos(t);

%% StateSpace Model
ss1 = ss(A, B, C, D);

%% Simulation
ZeroInputResponse = initial(ss1, xini, t);

[ImpulseInputResponse, t1] = impulse(ss1);
[StepInputResponse, te] = step(ss1);

LinearResponse = lsim(ss1, u, t, xini);

%% Plots
figure
subplot(221), plot(t, ZeroInputResponse)
xlabel('time [seconds]')
title('Zero Input Response')

subplot(222), plot(t1, ImpulseInputResponse)
xlabel('time [seconds]')
title('Impulse Input Response')

subplot(223), plot(te, StepInputResponse)
xlabel('time [seconds]')
title('Step Input Response')

subplot(224), plot(t, LinearResponse)
xlabel('time [seconds]')
title('Cosine Input Response')