% Copyright 2025, All Rights Reserved
% Code by Lloyd Fung
%    For package ODR-BINDy

clear all, close all, clc
figpath = './figs/';
addpath(genpath('./'));

% Random Seed
rng(12);
% Set highest polynomial order of the combinations of polynomials of the state vector
% polyorder = 3;
polyorder = 2;
% Set 1 to set yout = [yout sin(k*yin) cos(k*yin)];
usesine = 0; % For old poolData only
% Set the parameters of the Lorentz system
sigma = 10;  % Lorenz's parameters (chaotic)
beta = 8/3;
rho = 28;

% Set the number of columns of x (dimension of system)
D = 3;

%% Display ground truth
if polyorder == 2
    Xi_truth = zeros(10,3);
elseif polyorder == 3
    Xi_truth = zeros(20,3);
end    
% dy = [ sigma*(y(2)-y(1)) ; y(1)*(rho-y(3))-y(2) ; y(1)*y(2)-beta*y(3) ]
%               [   xdot , ydot , zdot  ]
Xi_truth(2,:) = [ -sigma ,  rho ,     0 ]; % y(1)
Xi_truth(3,:) = [ +sigma ,   -1 ,     0 ]; % y(2)
Xi_truth(4,:) = [      0 ,    0 , -beta ]; % y(3)
Xi_truth(6,:) = [      0 ,    0 ,     1 ]; % y(1) * y(2)
Xi_truth(7,:) = [      0 ,   -1 ,     0 ]; % y(1) * y(3)
disp('Ground Truth');
poolDataLIST({'x','y','z'},Xi_truth,D,polyorder,usesine);

%% generate Data
x0= [-8,8,27];%[-1,6,15]; %=[2.6,5,8];  % Initial condition
% Integrate in time using ode89
t_final=5;
dt=0.01;
tspan=dt:dt:t_final; % BINDy is best used in low data (can be lower than # of candidate func)
N = length(tspan);

ODEoptions = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,D));
[t,x_clean]=ode89(@(t,x) lorenz(t,x,sigma,beta,rho),tspan,x0,ODEoptions);

%% Compute Derivative and Add noise
%  Set size of random noise
eps_x = std(x_clean(:))*0.2; %Linear method work as well when noise is >2.0

% Actual dx from clean data. For reference only
dx_clean = NaN(length(x_clean),D);
% Calculate the derivative dx at each point along the trajectory
for i=1:length(x_clean)
    dx_clean(i,:) = lorenz(0,x_clean(i,:),sigma,beta,rho);
end

% Generating noisy data
x = x_clean + eps_x*randn(size(x_clean));

%% ODR-BINDy Greedy
% Library 
if polyorder == 2
    Libs.M = 10;
    Libs.Theta_fun = @(X)Polynomial3D2O(X);
    Libs.dTheta_fun = @(X)Polynomial3D2Od(X);
    Libs.ddTheta_fun = @(X)Polynomial3D2Odd(X);
    Libs.dddTheta_fun_f = @(X,p,mask)Polynomial3D2Oddd_f(X,p,mask);
elseif polyorder == 3
    Libs.M = 20;
    Libs.Theta_fun = @(X)Polynomial3D3O(X);
    Libs.dTheta_fun = @(X)Polynomial3D3Od(X);
    Libs.ddTheta_fun = @(X)Polynomial3D3Odd(X);
    Libs.dddTheta_fun_f = @(X,p,mask)Polynomial3D3Oddd_f(X,p,mask);
end

% FD object
ODR_int_pt=4;
[IMat,DMat]=FD(N,ODR_int_pt,dt,false);
TimeDiffObj.IMat = IMat;
TimeDiffObj.DMat = DMat;
TimeDiffObj.t = t;

% Hyperparam: Noise standard deviation and prior standard deviation
HyperObj.SigmaX = eps_x*ones(size(x));
HyperObj.SigmaY = 1e-4*ones(size(IMat,1),D);
HyperObj.SigmaP = 1e2*ones(Libs.M,3);

% Options (for plotting, output level, etc.)
ODRopts.PlotXout = true;
ODRopts.VerboseLevel = 2;
ODRopts.SaveProgress = false;
ODRopts.LoadProgressFile = "";% "ODRProgress.mat";

tic
[Xi_ODR,J_Evi_out,X_denoise]=ODR_BINDy_Greedy(x,Libs,TimeDiffObj,HyperObj,ODRopts);
% [Xi_ODR,J_Evi_out,X_denoise]=ODR_BINDy_Greedy_Parallel(x,Libs,TimeDiffObj,HyperObj,ODRopts);
disp('From Bayesian Regression (Greedy+ ODR)');
poolDataLIST({'x','y','z'},Xi_ODR,D,polyorder,usesine);
toc

return

%% LORENZ for T
savefig=true;
tA = tspan;
tO = tspan;
xA = x_clean;
xO = X_denoise;

%% Noise and DeNoise Visualisation
hf=figure('Position',[100 100 1250 450]);
a0=subplot(3,6,[1 2 7 8 13 14]);

hold on;
plot3(x(:,1),x(:,2),x(:,3),'rx:',LineWidth=0.5);
plot3(x_clean(:,1),x_clean(:,2),x_clean(:,3),'k-',LineWidth=2.0);
plot3(X_denoise(:,1),X_denoise(:,2),X_denoise(:,3),'--',...
    LineWidth=1.5);
view(27,16);
grid on;
set(gca,'TickDir','none','FontSize',12);
% set(get(gca, 'XAxis'), 'Visible', 'off');
% set(get(gca, 'YAxis'), 'Visible', 'off');
% set(get(gca, 'ZAxis'), 'Visible', 'off');
xlabel('$$x$$','Interpreter','latex','FontSize',14)
ylabel('$$y$$','Interpreter','latex','FontSize',14)
zlabel('$$z$$','Interpreter','latex','FontSize',14)
axis([-25 25 -35 35 0 60])
legend('Data','Truth','ODR-BINDy','Location','northeast','FontSize',12);
a0.Position=[0.04,0.130833333333333,0.3,0.794166666666666];

%% Lorenz, dynamo view
a1=subplot(3,6,[3 4 5 6]); hold on;
plot(t,x(:,1),'rx');
plot(tA,xA(:,1),'k-','LineWidth',2);
plot(tO,xO(:,1),'--','LineWidth',1.5)
grid on;
xlim([tspan(1) tspan(end)])
xticklabels({})
set(gca,'FontSize',12)
% xlabel('$$t$$','Interpreter','latex','FontSize',18)
ylabel('$$x$$','Interpreter','latex','FontSize',14)
a1.Position = [0.4,0.65,0.55,0.25];

a2=subplot(3,6,[9 10 11 12]); hold on;
plot(t,x(:,2),'rx');
plot(tA,xA(:,2),'k-','LineWidth',2);
plot(tO,xO(:,2),'--','LineWidth',1.5)
grid on;
xlim([tspan(1) tspan(end)])
xticklabels({})
set(gca,'FontSize',12)
% xlabel('$$t$$','Interpreter','latex','FontSize',18)
ylabel('$$y$$','Interpreter','latex','FontSize',14)
a2.Position=[0.4,0.38,0.55,0.25];

a3=subplot(3,6,[15 16 17 18]); hold on;
plot(t,x(:,3),'rx');
plot(tA,xA(:,3),'k-','LineWidth',2);
plot(tO,xO(:,3),'--','LineWidth',1.5)
grid on;
xlim([tspan(1) tspan(end)])

set(gca,'FontSize',12)
xlabel('$$t$$','Interpreter','latex','FontSize',14)
ylabel('$$z$$','Interpreter','latex','FontSize',14)
a3.Position = [0.4,0.11,0.55,0.25];
a3.YTick=0:20:60;
if savefig
    saveas(hf,[figpath 'Lorenz_dyn_eps' num2str(eps_x) '_dt' num2str(dt) '_tf' num2str(t_final) '.svg']);
    saveas(hf,[figpath 'Lorenz_dyn_eps' num2str(eps_x) '_dt' num2str(dt) '_tf' num2str(t_final) '.eps']);
    saveas(hf,[figpath 'Lorenz_dyn_eps' num2str(eps_x) '_dt' num2str(dt) '_tf' num2str(t_final) '.fig']);
end
