Qy = 1000;
Q = C' * Qy * C;
R = 1;

% LQR gain
K = lqr(A, B, Q, R);
figure;
pzmap(ss(A, B, C, 0)); hold on;
pzmap(ss(A - B*K, B, C, 0));
legend('Open-loop', 'LQR closed-loop');


% sys = ss(A, B, C, 0);
% rlocus(sys);
