function [Xi,J_Evi_out,X_denoised]=ODR_BINDy_Greedy_Parallel(Xdata,Libs,TimeDiffObj,HyperObj,options)
    %% ODR-BINDy - (Gaussian) Bayesian Inference and Model Selection by evidence maximisation with denoising via Orthogonal Distance Regression
    %
    % INPUT ARGUMENTS:
    % 
    %	Libs   Library of Candidate functions
    %
    %	Xdata		(orginal) Time-Series Data
    %
    %	TimeDiffObj.IMat		Matrix to compute the collocation point x
    % 
    %	TimeDiffObj.DMat		Matrix to compute the collocation point dx
    % 
    %	HyperObj.SigmaX    	Nosie (STD) in x
    %
    %   HyperObj.SigmaY      Noise (STD) in dx 
    %               (usually due to stochasticity or FD)
    %
    %   HyperObj.SigmaP   Prior sqrt(Variance)
    %
    %   options     Options for fine tuning noise iteration and turning
    %               on/off diagonstic graphs
    %
    % OUTPUT ARGUMENTS:
    %
    %   Xi     Output of estimated parameter (MAP point)
    %
    %   J_Evi_out  Negative Log-Evidence of the selected model
    %
    %   X_denoised Denoised x as a result of ODR

    % Copyright 2025, Code by Lloyd Fung
    %
    % This file is part of the ODR-BINDy package.
    % See "LICENSE" and "README.md" in the package for details about
    % license and copyright of the ODR-BINDy package.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Setting defaults for options
    if nargin <= 4
        options.MaxNoiseIter = 2; % Maximum noise iteration 
        options.StepTolerance = 1e-12; % Threshold for relative change in w for stopping noise iteration 
    end
    if ~isfield(options,'MaxNoiseIter')
        options.MaxNoiseIter = 2; % Maximum noise iteration 
    end
    if ~isfield(options,'StepTolerance')
        options.StepTolerance = 1e-12; % Threshold for relative change in w for stopping noise iteration 
    end
    if ~isfield(options,'LinUseDenoise')
        options.LinUseDenoise = true; % Threshold for relative change in w for stopping noise iteration 
    end
    if ~isfield(options,'PlotXout')
        options.PlotXout = true;
    end
    if ~isfield(options,'StartPtNum')
        options.StartPtNum = feature('numcores');
    end
    if ~isfield(options,'ContRunStop')
        options.ContRunStop = 2;
    end
    if ~isfield(options,'MaxFailedRun')
        options.MaxFailedRun = 5;
    end
    if ~isfield(options,'MaxInitialTrial')
        options.MaxInitialTrial = ceil(32/options.StartPtNum);
    end
    if ~isfield(options,'IniGuessUsePrev')
        options.IniGuessUsePrev = true; % Is it better to default this to false? seemed true is better??
    end
    if ~isfield(options,'UseGaussNewtonEst')
        options.UseGaussNewtonEst = true;
    end
    if ~isfield(options,'VerboseLevel')
        % 0: Disable all msg;
        % 1: At every removal of term
        % 2: Detailed
        % 3: All (including lsqnonlin internals)
        options.VerboseLevel = 2; 
    end
    if ~isfield(options,'SaveProgress')
        options.SaveProgress = false;
    elseif options.SaveProgress
        if ~isfield(options,'SaveProgressFileName')
            options.SaveProgressFileName = "ODRProgress";
        end
    end
    if ~isfield(options,'LoadProgressFile')
        options.LoadProgressFile="";
    end
    clear visualise_denoise;

    %% Initialisation
    % Extract the number of data points
    Nx = size(Xdata,1);
    Neq = size(TimeDiffObj.IMat,1);
    % Extract the number of terms in the Theta library
    M = Libs.M;
    % Extract the number of state variables
    ND = size(Xdata,2);
    
    % Boolean Mask
    mask = true(M*ND,1);
    NpAll = M*ND;

    Theta_fun = Libs.Theta_fun;
    dTheta_fun = Libs.dTheta_fun;

    Data = struct('Neq',Neq,'Nx',Nx,'ND',ND,'M',M,'Xdata',Xdata,'SigmaBeta',HyperObj.SigmaP,...
        'Theta_fun',Theta_fun,'dTheta_fun',dTheta_fun,'IMat',TimeDiffObj.IMat,'DMat',TimeDiffObj.DMat,...
        'SigmaX',HyperObj.SigmaX,'SigmaY',HyperObj.SigmaY);
    % Initialisation lsqnonlin options
    lsqnonlinopt = optimoptions('lsqnonlin','MaxFunctionEvaluations',1e6,'MaxIterations',1e2, ...
        'Algorithm','trust-region-reflective','FunctionTolerance',5e-8,'StepTolerance',options.StepTolerance,...
        'Display','off','SpecifyObjectiveGradient',true);
    lsqnonlinopt2 = optimoptions('lsqnonlin','MaxFunctionEvaluations',1e6,'MaxIterations',5e2, ...
        'Algorithm','trust-region-reflective','FunctionTolerance',5e-8,'StepTolerance',options.StepTolerance,...
        'Display','off','SpecifyObjectiveGradient',true);

    if options.VerboseLevel >=3
        lsqnonlinopt.Display = 'iter';
        lsqnonlinopt2.Display = 'iter';
    end
    % Initialise Initial Guesser with ensembling
    InitGuessOpt = struct('EnsembleSampleSize',100,'DataSelectRatio',1.0,'Bragging',true);
    InitGuessOpt2 = struct('EnsembleSampleSize',500,'DataSelectRatio',1.0,'Bragging',true);

    if options.PlotXout
        visualise_denoise(Xdata,Xdata);
    end
    %% Loop through each dimension
    % Initialise temporary array for storage
        w_j=cell(NpAll-ND,1);
     Xout_j=cell(NpAll-ND,1);
    J_Evi_j=Inf(1,NpAll-ND);

    % Greedy search: iterate by removing terms one-by-one in the model
    j=1;
    ContRun_flag=0;  
    proceed_flag=false;
    Failed_run = 0;
    J = Inf;
    % Initial Trial: Get initial guess of Xout to speed up iterations 
    if (options.LoadProgressFile=="") 
        while  J==Inf && Failed_run < options.MaxInitialTrial 
            [w,J,Xout_ini,exitflag_loc]=ODR_BINDy_Regression_MultiStart(Data,Libs,mask,lsqnonlinopt2,InitGuessOpt2,options);
    
            Failed_run = Failed_run +1;
            lsqnonlinopt2.MaxIterations=lsqnonlinopt2.MaxIterations*2;
            if options.VerboseLevel >=2
                disp(w);
            end
        end
    else
        exitflag_loc = 1;
    end

    if exitflag_loc <=0
        error('Initial Trial failed');
    else
        Failed_run = 0;
    end
            
    % Visualise
    if options.PlotXout && (options.LoadProgressFile=="") 
        visualise_denoise(Xout_ini);
    end
    if options.VerboseLevel >=1 && (options.LoadProgressFile=="") 
        disp('Initial guess of parameters');
        disp(w);
    end

    if ~(options.LoadProgressFile=="")
        StartPtNum = options.StartPtNum;
        load(options.LoadProgressFile); % Will overwrite j ContRun_flag Failed_run
        options.PlotXout = false; % Have to disable plotting. Figure handle is lost.
        Neq=Data.Neq;
        Nx=Data.Nx;
        ND=Data.ND;
        M=Data.M;
        NpAll=ND*M;
        options.StartPtNum=StartPtNum;
    end
    %% Greedy search: iterate by removing terms one-by-one in the model
    while j<=(NpAll-ND) && ContRun_flag < options.ContRunStop && Failed_run < options.MaxFailedRun
        Np = NpAll-j;
        J_Evi_i=Inf(1,Np+1);
        Xout_i=zeros(Nx,ND,Np+1);
        w_i=zeros(M,ND,Np+1);
        
        % Loop through all terms one can remove from current model
        m = (1:NpAll)';
        m = m(mask);
        parfor (ii = 1:Np+1, options.StartPtNum)
            if options.VerboseLevel >=2
                disp([num2str(ii) 'th variable zeroed | ' num2str(Np) ' variables total']);
            end
            i = m(ii);
            % Remove one term (the ith term from the current list of
            % indices m)
            mask_loc = mask;
            mask_loc(i) = false;
            if all(any(reshape(mask_loc,M,ND),1)) % Checking to make sure at least one term remains

                [w,J_Evi,Xout]=ODR_BINDy_Regression(Data,Libs,mask_loc,...
                    lsqnonlinopt,InitGuessOpt,options,Xout_ini);

                J_Evi_i(ii) = J_Evi;
                w_i(:,:,ii) = w;
                Xout_i(:,:,ii) = Xout;
            else
                J_Evi_i(ii) = Inf;
            end
        end

        if all(J_Evi_i==Inf)
            mask_= mask; % Run has failed. Re-run the previous initial trial to get better Xout_ini
            Failed_run = Failed_run + 1;
        else
            % Select the index with the MINIMUM J (maximum posterior likelihood).
            % This means: "when you drop the polynomial term with this index, the remaining model
            % has higher evidence than all the other alternative models."
            [~,i_sorted]=sort(J_Evi_i);
    
            % Remove this index from the indices that can be selected next
            % (avoid making no active term in each column)
            iii=1;
            mask_next = false(M,ND);
            while ~all(any(mask_next,1))
                mask_next = reshape(mask,M,ND);
                mask_next(m(i_sorted(iii))) = false;
                iii=iii+1;
            end
            mask_ = reshape(mask_next,M*ND,1);
            iii = iii-1;
            i_min = i_sorted(iii);
            % disp(w_i(:,:,i_min));

            % Reset Initial Trial Call
            proceed_flag=true;
            Failed_run = 0;
        end

        % Initial Trial: Get initial guess of Xout to speed up iterations
        if options.IniGuessUsePrev && proceed_flag
            [w_ini,J_ini,Xout_ini,exitflag_loc]=ODR_BINDy_Regression_MultiStart(Data,Libs,mask_,lsqnonlinopt2,InitGuessOpt2,options,Xout_i(:,:,i_min));
        else
            [w_ini,J_ini,Xout_ini,exitflag_loc]=ODR_BINDy_Regression_MultiStart(Data,Libs,mask_,lsqnonlinopt2,InitGuessOpt2,options);
        end
        if exitflag_loc >0 || J_ini==Inf
            Failed_run = Failed_run +1;
        end
        proceed_flag = proceed_flag && (exitflag_loc >0);

        % Visualise
        if options.PlotXout
            visualise_denoise(Xout_ini);
        end
        if options.VerboseLevel >=1
            disp(['After removal of the ' num2str(j) 'th term']);
            disp(w_j{j});
        end

        if proceed_flag
            if J_ini < J_Evi_i(i_min)
                % Storing values to array
                 Xout_j{j} = Xout_ini;
                    w_j{j} = w_ini; % History of parameters of terms selected
                J_Evi_j(j) = J_ini; % History of resulting neg. log-evidence
                
            else
                % Storing values to array
                 Xout_j{j} = Xout_i(:,:,i_min);
                    w_j{j} = w_i(:,:,i_min); % History of parameters of terms selected
                J_Evi_j(j) = J_Evi_i(i_min); % History of resulting neg. log-evidence
                Xout_ini = Xout_i(:,:,i_min);
            end
            disp(w_j{j});
            if j > 2
                if J_Evi_j(j)>=J_Evi_j(j-1)
                    ContRun_flag=ContRun_flag+1; % Raise the counter if J_evi start rising
                else
                    ContRun_flag=0; % Reset counter if J_evi drops
                end
            end

            mask=mask_;
            j=j+1;
            proceed_flag = false; % Reset
        end
        if options.SaveProgress
            save(options.SaveProgressFileName,...
                "options","TimeDiffObj","HyperObj","Data",...
                "j","mask","Xout_ini",...
                "ContRun_flag","Failed_run","proceed_flag",...
                "Xout_j","w_j","J_Evi_j",...
                '-v7.3');
        end
    end
    % Selects the MODEL with the minimum J, which is the maximum evidence 
    [J_Evi_out,j_min]=min(J_Evi_j);
    % Format the parameter estimate back into a sparse matrix
    Xi = w_j{j_min};
    X_denoised = Xout_j{j_min};

    if options.PlotXout
        visualise_denoise(X_denoised)
    end
end

function visualise_denoise(Xpred,varargin)
    persistent ax ND
    if nargin > 1
        Xdata = varargin{1};
    end
    if isempty(ax) || nargin > 1
        % Initialise 
        ND = size(Xdata,2);

        figure;
        clr=colororder("gem");
        
        if ND <= size(clr,1)
            colororder(clr(1:ND,:));
        end

        plot(Xdata,'x');hold on;

        ax = plot(Xpred);hold off;
        
        drawnow();
        set(gca,'YLimMode','manual');
    else
        % Update plot
        for plot_i = 1:ND
            ax(plot_i).YData=Xpred(:,plot_i);
        end
        drawnow();
    end
end