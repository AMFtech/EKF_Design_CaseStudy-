%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate 100 Sensory Observation Data Set, each contains 500 points in
% polar coordinate system, given raw trajectory in Cartesian Coordinate
% SohwPlot:         is a flag to show a plot of 
% SonarTrjcDataSet: m*n*100 Matrix, each m*n Matrix contains a single...
%                   trajectory: [RHO_k, THETA_k] 
%% 
function [SonarDataSet, SonarDataSetCart] = SonarDataSetGen(ShowPlot)

global Sig_rho;
global Sig_theta;

if(nargin == 0)
    ShowPlot = 0;
end

load object_trajectory.mat
RealPos = x_traj_pos;

[m, n] = size(RealPos);
X_RealPos = RealPos(1,:);
Y_RealPos = RealPos(2,:);

SonarDataSet        = zeros(m,n,100);
SonarDataSetCart    = zeros(m,n,100);

for k = 1:100
    [THETA, RHO] = cart2pol(X_RealPos, Y_RealPos);
    THETAz  = THETA + randn(1, n) * Sig_theta;  % SONAR RHO
    RHOz    = RHO   + randn(1, n) * Sig_rho;    % SONAR Theta
    [Zx, Zy] = pol2cart(THETAz, RHOz);
    
    SonarDataSet(1,:,k) = RHOz;
    SonarDataSet(2,:,k) = THETAz;
    
    if(nargout == 2)
        SonarDataSetCart(1,:,k) = Zx;
        SonarDataSetCart(2,:,k) = Zy;
    end    
    
    if (ShowPlot == 1)
        plot(Zx,Zy, '--b', 'LineWidth', 1);
        hold on;
    end
end

if (ShowPlot == 1)
    plot(X_RealPos, Y_RealPos, 'r', 'LineWidth',3 );
    title('Red: Real Trajecory, Blue: Noisy Trajectories for all data set')
    xlabel('X position in Cartesian Coordinate', 'FontSize', 12, ...
        'FontWeight','bold');
    ylabel('Y position in Cartesian Coordinate', 'FontSize', 12, ...
        'FontWeight','bold');    
end
