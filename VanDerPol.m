% Copyright 2025, All Rights Reserved
% Code by Lloyd Fung
%    For package ODR-BINDy

clear all, close all, clc
figpath = './figs/';
addpath(genpath('./'));

% Random Seed
rng(129);
% Set highest polynomial order of the combinations of polynomials of the state vector
polyorder = 3;
% Set 1 to set yout = [yout sin(k*yin) cos(k*yin)];
usesine = false; % For old poolData only
% Set the parameters of the Van Der Pol system
param.beta  = 2.5;

% Set the number of columns of x (dimension of system)
D = 2;

%% Display ground truth
Xi_truth = zeros(10,2);
Xi_truth(3,1)=1;
Xi_truth(2,2)=-1;
Xi_truth(3,2)=param.beta ;
Xi_truth(8,2)=-param.beta ;
disp('Ground Truth');
poolDataLIST({'x','y'},Xi_truth,D,polyorder,usesine);

%% Generate Data
x0=[-2 1];
% Integrate in time using ode89
t_final=10;
dt=0.05; %0.025
tspan=dt:dt:t_final; % BINDy is best used in low data (can be lower than # of candidate func)
N = length(tspan);

ODEoptions = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,2),InitialStep=1e-6);
[t,x_clean]=ode89(@(t,x) vanderpol(t,x,param),tspan,x0,ODEoptions);

%% Compute Derivative and Add noise
%  Set size of random noise
eps_x = mean(rms(x_clean))*0.2;

% Actual dx from clean data. For reference only
dx_clean = NaN(length(x_clean),D);
% Calculate the derivative dx at each point along the trajectory
for i=1:length(x_clean)
    dx_clean(i,:) = vanderpol(t(i),x_clean(i,:),param);
end

% Generating noisy data
x = x_clean + eps_x*randn(size(x_clean));

%% ODR-BINDy Greedy
% Library 
Libs.M = 10;
Libs.Theta_fun = @(X)Polynomial2D3O(X);
Libs.dTheta_fun = @(X)Polynomial2D3Od(X);
Libs.ddTheta_fun = @(X)Polynomial2D3Odd(X);
Libs.dddTheta_fun_f = @(X,p,mask)Polynomial2D3Oddd_f(X,p,mask);

% FD object
ODR_int_pt = 4;
[IMat,DMat]=FD(N,ODR_int_pt,dt,false);
TimeDiffObj = struct('t',tspan);
TimeDiffObj.IMat = IMat;
TimeDiffObj.DMat = DMat;

% Hyperparam: Noise standard deviation and prior standard deviation
HyperObj = struct('SigmaX',eps_x*ones(size(x)));
HyperObj.SigmaY = 1e-4*ones(size(IMat,1),D); 
HyperObj.SigmaP = 1e1*ones(Libs.M,3);

ODRopts.PlotXout = true;
ODRopts.VerboseLevel = 2;
tic
% [Xi_ODR,~,X_denoise]=ODR_BINDy_Greedy(x,Libs,TimeDiffObj,HyperObj,ODRopts);
[Xi_ODR,J_Evi_out,X_denoise]=ODR_BINDy_Greedy_Parallel(x,Libs,TimeDiffObj,HyperObj,ODRopts);
disp('From Bayesian Regression (Greedy+ ODR)');
poolDataLIST({'x','y'},Xi_ODR,D,polyorder,usesine);
toc
