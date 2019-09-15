%% The input is the Polar Coordinate sensory read out from Sonar Sensor Zk
%   Inputs: X_EKF:          Estimated EKF state vectors     
%           K_Den_Factor:   Denominator factor of Kalman Gain
%           
% ASSUMPTION:   THE 1st ROW OF Zk IS RHO, THE 2nd IS THETA

function [X_EKF, K_Den_Factor, Kk_NORM ] = myEKF(Zk_Sonar)

global q;           % Process Covariance parameter
global Sig_rho;     % RHO variance in meter
global Sig_theta;   % Theta variance in rad

r = 1;

[m, n] = size(Zk_Sonar);
%% Matrices Initialization 

% Covariance Matrix of process noise... (Velocity Noise)
Ts = 0.1;           % sampling time, dt = 0.1 sec
Q = q * [ Ts^3/3, 0,      Ts^2/2, 0;
          0,      Ts^3/3, 0,      Ts^2/2; 
          Ts^2/2, 0       Ts,     0;
          0,      Ts^2/2, 0       Ts];

% State Transition Matrix
F = [ 1, 0, Ts, 0;
      0, 1, 0   Ts; 
      0, 0, 1,  0
      0, 0, 0,  1];

% Initial State of State Covariance Matrix
Pk_1 = [ 0, 0, 0, 0;
         0, 0, 0, 0;
         0, 0, 0, 0;
         0, 0, 0, 0;];
    
% Observation Covariance Matrix
R = r*[ Sig_rho^2, 0;  
      0,         Sig_theta^2];

% Initial State Vector
Xk_1 = [50; 50; 0; 0];

% Stored State Estimate, outcome of EKF
X_EKF = zeros(4, n);
K_Den_Factor = zeros(1, n);
Kk_NORM = zeros(1, n);

% Recursive Kalman Operation-----------------------------------------------
for i = 1 : 500
    % Prediction phase --------------------------------
    Xp = F*Xk_1;            % Initial state estimate
    Pp = F*Pk_1*F' + Q;     % Initial Process Covariance Matrix 
    %---------------------------------------------------
 
    % Measurement update ------------------------------
    [THETAp, RHOp] = cart2pol(Xp(1), Xp(2));
    
    % Compute the Jacobian of Zk = h(Xk;Vk);
    Hk = [cos(THETAp),       sin(THETAp),      0, 0;
          -sin(THETAp)/RHOp, cos(THETAp)/RHOp, 0, 0];
    K = (Pp*Hk')/(Hk*Pp*Hk' + R);      % Kalman Gain
    dZk = [Zk_Sonar(1,i); Zk_Sonar(2,i)] - [RHOp ; THETAp] ;
 
    Xk = Xp + K * dZk;
    Pk = Pp - K * Hk * Pp;
    X_EKF(:,i) = Xk;
    %---------------------------------------------------
    
    % Update the initial values for next state estimation
    Xk_1 = Xk;
    Pk_1 = Pk;
    
    % This factor reflects the relative weight of Process Noise vs Sensory
    % noise: Greater the factor, Stronger the process noise.
    K_Den_Factor(i) = norm( Hk*Pp*Hk') / norm(R);
%     K_Den_Factor(i) = norm( eye(4) - K*Hk ) / norm(K);
    Kk_NORM(i) = norm(Pk);

end

end
