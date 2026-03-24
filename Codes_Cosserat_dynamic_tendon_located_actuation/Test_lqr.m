
% LQR setup: output-weighted cost
Qy = 1;
Q = C' * Qy * C;
R = 1;

% LQR gain
K = lqr(A, B, Q, R);

% Feedforward gain for reference tracking
%N = inv(C * inv(A - B*K) * B);
N = -inv(C * inv(A - B*K) * B);  % Remove minus sign if it's there

% Reference output
y_ref = 0.01;

% Time and initial condition
t_span = [0 10];
x0 = zeros(6,1);

% Dynamics function for ode45 (with clipped input)
lqr_dynamics = @(t, x) A * x + B * max(0, -K * x + N * y_ref);  
% u = max(0, -Kx + Ny_ref)
% Simulate with ode45
[t, x] = ode45(lqr_dynamics, t_span, x0);

% Compute output and control signal
y = x * C';
u = max(0, -x * K' + N * y_ref);  % Control input (tension), clipped
%C*(A - BK)^(-1)*B
% Plot output tracking
figure;
subplot(2,1,1);
plot(t, y, 'b', 'LineWidth', 2);
yline(y_ref, 'r--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Output y(t)');
title('LQR Output Tracking (with Positive Cable Tension)');
legend('y(t)', 'y_{ref}');
grid on;

% Plot control signal
subplot(2,1,2);
plot(t, u, 'k', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Control Input u(t)');
title('Cable Tension (Clipped at u ≥ 0)');
grid on;
