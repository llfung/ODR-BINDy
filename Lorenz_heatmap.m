% Copyright 2025, All Rights Reserved
% Code by Lloyd Fung
%    For package ODR-BINDy
%
% Run simulations over noise level and data length for heatmap plot
%

addpath(genpath('./'));

%% sweep over a set of noise levels and data length to generate heatmap plots
% noise level
epsL = 0.025:0.025:0.4;

% simulation time
tEndL = 4.0:1.0:10;

% at each noise level and simulation time, nTest different instantiations of noise are run (model errors and success rate are then averaged for plotting)
nTest1 = 128; % generate models nTest1 times for SINDy

%% Settings
% Set highest polynomial order of the combinations of polynomials of the state vector
% polyorder = 3;
polyorder = 2;

% Set 1 to set yout = [yout sin(k*yin) cos(k*yin)]  (For poolDataLIST)
usesine = 0;

% Set the number of dimensions
D = 3;

%% hyperparameters
PparamV_ODR= 100^2;% Arbitrary large variance with zero mean for all coefficients
SigmaY_ODR = 1e-3;
ODR_int_pt = 6;

%% Build library of nonlinear terms as functions
if polyorder == 2
    libs.M = 10;
    libs.Theta_fun = @(X)Polynomial3D2O(X);
    libs.dTheta_fun = @(X)Polynomial3D2Od(X);
    libs.ddTheta_fun = @(X)Polynomial3D2Odd(X);
    libs.dddTheta_fun_f = @(X,p,mask)Polynomial3D2Oddd_f(X,p,mask);
elseif polyorder == 3
    libs.M = 20;
    libs.Theta_fun = @(X)Polynomial3D3O(X);
    libs.dTheta_fun = @(X)Polynomial3D3Od(X);
    libs.ddTheta_fun = @(X)Polynomial3D3Odd(X);
    libs.dddTheta_fun_f = @(X,p,mask)Polynomial3D3Oddd_f(X,p,mask);
end

%% common parameters, true Lorenz system, signal power for noise calculation

% generate synthetic Lorenz system data
sigma = 10;  % Lorenz's parameters (chaotic)
beta = 8/3;
rho = 28;
x0 = [-8 8 27]';


% set common params
tol_ode = 1e-10;         % set tolerance (abs and rel) of ode89
options = odeset('RelTol',tol_ode,'AbsTol',tol_ode*ones(1,length(x0)));

% time step
dt = 0.01;

% get true Lorenz system for comparison
Xi_truth = zeros(libs.M ,D);
Xi_truth( 2,:) = [-sigma rho 0];
Xi_truth( 3,:) = [ sigma -1  0];
Xi_truth( 4,:) = [ 0      0  -beta];
Xi_truth( 6,:) = [ 0      0  1];
Xi_truth( 7,:) = [ 0     -1  0];

% signal power for noise calculation
[~,x_]=ode89(@(t,x) lorenz(t,x,sigma,beta,rho),dt:dt:10,x0,options);

signal_power = std(x_(:));

%% general parameters
saveTrue = 0;

%% Run the Loop
% Initialisation
nWrongTermsODR = zeros(length(epsL),length(tEndL),nTest1);
modelErrorODR = zeros(length(epsL),length(tEndL),nTest1);
successODR = zeros(length(epsL),length(tEndL),nTest1);

% Loop
for ieps = 1:length(epsL)
    for idt = 1:length(tEndL)
        %% Iterate through each noise and data length
        noise_ratio= epsL(ieps);
        tEnd = tEndL(idt);

        tspan = dt:dt:tEnd;

        [t,x_clean]=ode89(@(t,x) lorenz(t,x,sigma,beta,rho),tspan,x0,options);

        dx_clean = zeros(size(x_clean));
        for i=1:size(x_clean,1)
            dx_clean(i,:) = lorenz(t(i),x_clean(i,:),sigma,beta,rho);
        end

            nWrongTermsODR_temp = zeros(nTest1,1);
            modelErrorODR_temp  = zeros(nTest1,1);
            successODR_temp = zeros(nTest1,1);
            VectorErrorODR_temp = zeros(nTest1,1);
            NoiseIDErrorODR_temp = zeros(nTest1,1);

        parfor ii = 1:nTest1
            % set rnd number for randomness
            rng("default");

            % add noise
            eps_x = noise_ratio*signal_power;
            noise = normrnd(0,eps_x,size(x_clean));
            x = x_clean + noise;
            Theta_clean = poolData(x_clean,D,polyorder,usesine);

                %% ODR-BINDy Greedy Settings
                N=length(tspan);
                % FD object
                [IMat,DMat]=FD(N,ODR_int_pt,dt,false);
                TimeDiffObj = struct('t',tspan);
                TimeDiffObj.IMat = IMat;
                TimeDiffObj.DMat = DMat;

                % Hyperparam: Noise standard deviation and prior standard deviation
		        HyperObj = struct('SigmaX',eps_x*ones(size(x)));
                HyperObj.SigmaY = SigmaY_ODR*ones(size(IMat,1),D); 
                HyperObj.SigmaP = sqrt(PparamV_ODR)*ones(libs.M,3);


                %% ODR-BINDy
                ODRopts = struct('PlotXout',false,'VerboseLevel',0);
                [Xi_ODR,~,X_ODR]=ODR_BINDy_Greedy(x,libs,TimeDiffObj,HyperObj,ODRopts);

                %% store outputs
                if isempty(Xi_ODR)
                    nWrongTermsODR_temp(ii) = NaN;
                    modelErrorODR_temp(ii) = NaN;
                    successODR_temp(ii) = 0;
                    VectorErrorODR_temp(ii) = NaN;
                    NoiseIDErrorODR_temp(ii) = NaN;
                else
                    nWrongTermsODR_temp(ii) = sum(sum(abs((Xi_truth~=0) - (Xi_ODR~=0))));
                    modelErrorODR_temp(ii) = norm(Xi_ODR(:)-Xi_truth(:),'fro')/norm(Xi_truth(:),'fro');
                    successODR_temp(ii) = norm((Xi_truth~=0) - (Xi_ODR~=0))==0;
                    VectorErrorODR_temp(ii) = sum((dx_clean-Theta_clean*Xi_ODR).^2,"all")/sum(dx_clean.^2,"all");
                    NoiseIDErrorODR_temp(ii) = sum(((X_ODR-x_clean)-noise).^2,"all");
                end
        end
        %% Store outputs (outer-loop)
        nWrongTermsODR(ieps,idt,:) = nWrongTermsODR_temp;
        modelErrorODR(ieps,idt,:) = modelErrorODR_temp;
        successODR(ieps,idt,:) = successODR_temp;
        if saveTrue
            iii=(ieps-1)*length(tEndL)+idt;
            save([num2str(iii) '_Lorenz_SuccessRate_heatmap.mat'],...
                "nWrongTermsODR_temp","modelErrorODR_temp",...
                "successODR_temp","VectorErrorODR_temp","NoiseIDErrorODR_temp");
        end
    end
end
clearvars pc libs ans
save('Full_Lorenz_SuccessRate_heatmap.mat','-v7.3');
return