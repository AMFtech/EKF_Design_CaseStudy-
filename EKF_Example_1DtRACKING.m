% KALMAN FILTER EXMPLE FOR TRACKING SINGLE CART 
% CONSTANT VELOCITY ASSUMPTION IS MADE
 
clc
close all;
clear all;


dt = 1;         % dt = 1 sec
P = [ 20, 0;
       0, 10; ];
X = [ 0;
      0; ];
R = [ 1; ];
F = [ 1, dt;
      0, 1; ];
H = [ 1, 0; ];
Q = (1e-6) * [ (dt^4)/4, (dt^3)/2;
                 (dt^3)/2, dt^2; ];
pos = 0;
nIter = 200;
velocities = zeros(1, nIter);
positions = zeros(1, nIter);
truePositions = zeros(1, nIter);
KalmanGain = zeros(2, nIter);

for i = 1 : nIter
    Xp = F*X;
    Pp = F*P*F' + Q;
 
    K = (Pp*H')/(H*Pp*H' + R);
    y = ( pos + randn()*5 ) - H*Xp;
 
    X = Xp + K * y;
    P = Pp - K * H * Pp;
 
    KalmanGain(:,i) = K;
    NormK(i) = norm(K);
    
    velocities(i) = X(2);
    positions(i) = X(1);
    truePositions(i) = pos;
 
    pos = pos + 0.1;
end
 
plot(linspace(0, nIter, nIter), velocities, 'r', 'LineWidth',2);
title('Estimated Velocity');
xlabel('s')
ylabel('m/s');
 
figure;
plot(linspace(0, nIter, nIter), positions, 'r', 'LineWidth', 2);
hold on;
plot(linspace(0, nIter, nIter), truePositions, '--b', 'LineWidth', 2);
title('Estimated Position');
xlabel('s')
ylabel('m');

figure;
plot(positions, velocities, '--r', 'MarkerSize', 5, 'LineWidth',1, ...
    'MarkerFaceColor', 'b');
title('Kalman Gain eterminant');
xlabel('KalmanGain for x')
ylabel('KalmanGain for dx');

figure;
plot(KalmanGain(1,:), KalmanGain(2,:), '--o', 'Color', 'red', 'MarkerSize', 5, ...
    'LineWidth',1, 'MarkerFaceColor', 'b');
title('Kalman Gain eterminant');
xlabel('KalmanGain for x')
ylabel('KalmanGain for dx');

