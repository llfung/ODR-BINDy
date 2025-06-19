% Copyright 2025, All Rights Reserved
% Code by Lloyd Fung
%    For package ODR-BINDy

clear all, close all, clc
figpath = './figs/';
addpath(genpath('./'));

% Random Seed
rng(129);
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
t_final=10; %4 is better than 40?
dt=t_final/1000;
tspan=dt:dt:t_final; % BINDy is best used in low data (can be lower than # of candidate func)
N = length(tspan);

ODEoptions = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,D));
[t,x_clean]=ode89(@(t,x) lorenz(t,x,sigma,beta,rho),tspan,x0,ODEoptions);

%% Compute Derivative and Add noise
%  Set size of random noise
eps_x = 16*0.2; %Linear method work as well when noise is >2.0

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

%% FIGURE 1:  LORENZ for T\in[0,20]
savefig=false;
tspan = [0 4];
[tA,xA]=ode89(@(t,x)lorenz(t,x,sigma,beta,rho),tspan,x0,ODEoptions);   % true model
[tO,xO]=ode89(@(t,x)sparseGalerkin(t,x,Xi_ODR,polyorder,usesine),tspan,x0,ODEoptions);  % SINDY approximate

%%
f1=figure('Position',[200,259,390,413]);
% subplot(1,3,1)
plot3(xA(:,1),xA(:,2),xA(:,3),'LineWidth',1.5);hold on;
plot3(x(:,1),x(:,2),x(:,3),'x');
view(27,16)
grid on
xlabel('x','FontSize',13)
ylabel('y','FontSize',13)
zlabel('z','FontSize',13)
set(gca,'FontSize',13)
% title('Orignal');
legend('Truth','Data','Location','northeast');
axis([-25 25 -35 35 0 60]);
if savefig
    saveas(f1,[figpath 'Lorenz_Theo_eps' num2str(eps_x) '_dt' num2str(dt) '_tf' num2str(t_final) '.svg']);
    saveas(f1,[figpath 'Lorenz_Theo_eps' num2str(eps_x) '_dt' num2str(dt) '_tf' num2str(t_final) '.fig']);
end

f4=figure('Position',[1247,259,390,413]);
plot3(xO(:,1),xO(:,2),xO(:,3),'LineWidth',1.5); hold on;
plot3(x(:,1),x(:,2),x(:,3),'x');
view(27,16)
grid on
xlabel('x','FontSize',13)
ylabel('y','FontSize',13)
zlabel('z','FontSize',13)
% title('SINDy');
legend('ODR-BINDy','Location','northeast');
set(gca,'FontSize',13)
axis([-25 25 -35 35 0 60])
if savefig
    saveas(f4,[figpath 'Lorenz_ODR-BINDy_eps' num2str(eps_x) '_dt' num2str(dt) '_tf' num2str(t_final) '.svg']);
    saveas(f4,[figpath 'Lorenz_ODR-BINDy_eps' num2str(eps_x) '_dt' num2str(dt) '_tf' num2str(t_final) '.fig']);
end
%%
% Lorenz for t=20, dynamo view
f5=figure('Position',[100 100 1250 250]);
subplot(1,12,[1 2 3]); hold on;
plot(tA,xA(:,1),'k-','LineWidth',2);
plot(t,x(:,1),'kx');
plot(tO,xO(:,1),'--','LineWidth',1.5)
grid on;
xlim(tspan)
set(gca,'FontSize',14)
xlabel('$$t$$','Interpreter','latex','FontSize',18)
ylabel('$$x_1$$','Interpreter','latex','FontSize',18)


subplot(1,12,[5 6 7]); hold on;
plot(tA,xA(:,2),'k-','LineWidth',2);
plot(t,x(:,2),'kx');
plot(tO,xO(:,2),'--','LineWidth',1.5)
grid on;
xlim(tspan)
set(gca,'FontSize',14)
xlabel('$$t$$','Interpreter','latex','FontSize',18)
ylabel('$$x_2$$','Interpreter','latex','FontSize',18)


subplot(1,12,[9 10 11]); hold on;
plot(tA,xA(:,3),'k-','LineWidth',2);
plot(t,x(:,3),'kx');
plot(tO,xO(:,3),'--','LineWidth',1.5)
grid on;
xlim(tspan)
set(gca,'FontSize',14)
xlabel('$$t$$','Interpreter','latex','FontSize',18)
ylabel('$$x_3$$','Interpreter','latex','FontSize',18)
legend('Truth','Data','ODR-BINDy','Position',[0.847629796839729 0.82 0.144894986231761 0.13782991202346]);
% legend('Truth','Data','Bayesian-SINDy','Location','eastoutside');

% % Create textbox
% annotation(f5,'textbox',...
%     [0.0690067720090291 0.924630498533724 0.0420428893905192 0.0762463343108505],...
%     'String',{'(a)'},...
%     'LineStyle','none',...
%     'FontSize',20,...
%     'FontName','Times New Roman',...
%     'FontAngle','italic');
% 
% % Create textbox
% annotation(f5,'textbox',...
%     [0.632212189616249 0.927563049853372 0.0414785553047404 0.0762463343108505],...
%     'String',{'(c)'},...
%     'LineStyle','none',...
%     'FontSize',20,...
%     'FontName','Times New Roman',...
%     'FontAngle','italic');
% 
% % Create textbox
% annotation(f5,'textbox',...
%     [0.350045146726861 0.927563049853372 0.0420428893905192 0.0762463343108505],...
%     'String',{'(b)'},...
%     'LineStyle','none',...
%     'FontSize',20,...
%     'FontName','Times New Roman',...
%     'FontAngle','italic');

if savefig
    saveas(f5,[figpath 'Lorenz_dyn_eps' num2str(eps_x) '_dt' num2str(dt) '_tf' num2str(t_final) '.svg']);
    saveas(f5,[figpath 'Lorenz_dyn_eps' num2str(eps_x) '_dt' num2str(dt) '_tf' num2str(t_final) '.fig']);
end
%% Noise and DeNoise Visualisation
hf=figure('Position',[10 10 250 350]);hold on;
plot3(x(:,1),x(:,2),x(:,3),'r--',LineWidth=0.2);
plot3(x_clean(:,1),x_clean(:,2),x_clean(:,3),'k-',LineWidth=3.0);
view([391.185546875 15.41015625]);
grid on;
set(gca,'TickDir','none')
set(get(gca, 'XAxis'), 'Visible', 'off');
set(get(gca, 'YAxis'), 'Visible', 'off');
set(get(gca, 'ZAxis'), 'Visible', 'off');
axis([-20 20 -22.8 22 5 45]);


hf2=figure('Position',[260 10 250 350]);hold on;
plot3(x_clean(:,1),x_clean(:,2),x_clean(:,3),'k-',LineWidth=3.0);
plot3(X_denoise(:,1),X_denoise(:,2),X_denoise(:,3),'--',...
    LineWidth=2.0,Color=[0.0745098039215686 0.623529411764706 1]);
view([391.185546875 15.41015625]);
grid on;
set(gca,'TickDir','none')
set(get(gca, 'XAxis'), 'Visible', 'off');
set(get(gca, 'YAxis'), 'Visible', 'off');
set(get(gca, 'ZAxis'), 'Visible', 'off');
axis([-20 20 -22.8 22 5 45]);