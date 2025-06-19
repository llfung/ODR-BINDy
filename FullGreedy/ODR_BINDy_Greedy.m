function [Xi,J_Evi_out,X_denoised]=ODR_BINDy_Greedy(Xdata,Libs,TimeDiffObj,HyperObj,options)
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
    if ~isfield(options,'MinInitialTrial')
        options.MinInitialTrial = 4;
    end
    if ~isfield(options,'MaxInitialTrial')
        options.MaxInitialTrial = 16;
    end
    if ~isfield(options,'MaxFailedRun')
        options.MaxFailedRun = 2;
    end
    if ~isfield(options,'ContRunStop')
        options.ContRunStop = 2;
    end
    if ~isfield(options,'IniGuessUsePrev')
        options.IniGuessUsePrev = true; % Seems making this true speed up without much compromise on accuracy
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
    lsqnonlinopt = optimoptions('lsqnonlin','MaxFunctionEvaluations',1e6,'MaxIterations',100, ...
        'Algorithm','trust-region-reflective','FunctionTolerance',5e-8,'StepTolerance',options.StepTolerance,...
        'Display','off','SpecifyObjectiveGradient',true);
    lsqnonlinopt2 = optimoptions('lsqnonlin','MaxFunctionEvaluations',1e6,'MaxIterations',1000, ...
        'Algorithm','trust-region-reflective','FunctionTolerance',5e-8,'StepTolerance',options.StepTolerance,...
        'Display','off','SpecifyObjectiveGradient',true);

    if options.VerboseLevel >=3
        lsqnonlinopt.Display = 'iter';
        lsqnonlinopt2.Display = 'iter';
    end
    % Initialise Initial Guesser with ensembling
    InitGuessOpt = struct('EnsembleSampleSize',100,'DataSelectRatio',1.0,'Bragging',true);
    InitGuessOpt2 = struct('EnsembleSampleSize',1000,'DataSelectRatio',1.0,'Bragging',true);

    if options.PlotXout
        visualise_denoise(Xdata,Xdata);
    end
    %% Loop through each dimension
    % Initialise temporary array for storage
        w_j=cell(NpAll-ND,1);
     Xout_j=cell(NpAll-ND,1);
    J_Evi_j=Inf(1,NpAll-ND);

    % Initial Trial: Get initial guess of Xout to speed up iterations
    exitflag_loc=-1;
    J=Inf;
    Ini_trial_call = 0;
    while (Ini_trial_call < options.MinInitialTrial || exitflag_loc<=0 || J==Inf) && Ini_trial_call < options.MaxInitialTrial && (options.LoadProgressFile=="") 
        [w_,J_,Xout_ini_,exitflag_loc]=ODR_BINDy_Regression(Data,Libs,mask,lsqnonlinopt2,InitGuessOpt2,options);

        Ini_trial_call = Ini_trial_call +1;
        lsqnonlinopt2.MaxIterations = lsqnonlinopt2.MaxIterations*2; % might just be not enough iterations
        if J_< J && exitflag_loc > 0
            J = J_;
            w=w_;
            Xout_ini = Xout_ini_;
        end
        if options.VerboseLevel >=2
            disp(w_);
        end
        % 
        if options.PlotXout  
            visualise_denoise(Xout_ini_);
        end
    end

    %%
    % Exit run if the initial trial has failed
    if Ini_trial_call >= options.MaxInitialTrial && (options.LoadProgressFile=="") 
        disp('Number of Initial Trial has exceed MaxInitialTrial. Run failed');
        % Export the failed run result nonetheless
        Xi = w_;
        J_Evi_out = Inf;
        X_denoised = Xout_ini_;

        if options.PlotXout
            visualise_denoise(X_denoised);
        end
    else
        if options.VerboseLevel >=1 && (options.LoadProgressFile=="")
            disp('Initial guess of parameters');
            disp(w);
        end
        %% Continue run if the initial trial is successful
        % Greedy search initialisation of flags and counts
        j=1;
        ContRun_flag=0;  
        Failed_run = 0;
        if ~(options.LoadProgressFile=="")
            load(options.LoadProgressFile); % Will overwrite j ContRun_flag Failed_run
            options.PlotXout = false; % Have to disable plotting. Figure handle is lost.
            Neq=Data.Neq;
            Nx=Data.Nx;
            ND=Data.ND;
            M=Data.M;
            NpAll=ND*M;
        end
        %% Greedy search: iterate by removing terms one-by-one in the model
        while j<=(NpAll-ND) && ContRun_flag < options.ContRunStop && Failed_run < options.MaxFailedRun
            Np = NpAll-j;
            J_Evi_i=Inf(1,Np+1);
            Xout_i=zeros(Nx,ND,Np+1);
            w_i=zeros(M,ND,Np+1);
            
            %% Loop through all terms one can remove from current model
            m = (1:NpAll)';
            m = m(mask);
            for ii = 1:Np+1
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
    
                    if options.PlotXout
                        visualise_denoise(Xout);
                    end
    
                else
                    J_Evi_i(ii) = Inf;
                end
            end
            %% Post-processing the Greedy search result: Decision on term to remove
            RemovalDecision_flag = false; % Indicate a removal decision has been made
            while ~RemovalDecision_flag
                % Try to find a term to remove. Return proceed_flag=flase if can't
                [proceed_flag,mask_,i_min]=RemovalDecision(mask,J_Evi_i,m,ND); 
    
                %% Check chosen term to remove can indeed be removed (see convergence)
                % and get initial guess of Xout after removal of terms to speed up iterations
                if options.IniGuessUsePrev && proceed_flag
                    [w_ini_,J_ini_,Xout_ini_,exitflag_loc]=ODR_BINDy_Regression(Data,Libs,mask_,lsqnonlinopt2,InitGuessOpt2,options,Xout_i(:,:,i_min));
                else
                    exitflag_loc=-1; J_ini_ = Inf;
                end
                % If the Xout did not converge, or if we can't find a suitable term
                % to remove, redo the inital run without previous guess
                while (exitflag_loc<=0 || J_ini_==Inf) && Failed_run < options.MaxFailedRun
                    [w_ini_,J_ini_,Xout_ini_,exitflag_loc]=ODR_BINDy_Regression(Data,Libs,mask_,lsqnonlinopt2,InitGuessOpt2,options);
                    Failed_run = Failed_run +1;
                end
                if proceed_flag  && (Failed_run >= options.MaxFailedRun) % If removing the selected term leads to non-convergence, try another one
                    J_Evi_i(i_min) = Inf;
                    Failed_run = 0; % Reset Failed Count 
                else
                    RemovalDecision_flag = true; % Decision is finalised to proceed.
                end
            end
                
            if proceed_flag
                if J_ini_ < J_Evi_i(i_min)
                    % Storing values to array
                     Xout_j{j} = Xout_ini_;
                        w_j{j} = w_ini_; % History of parameters of terms selected
                    J_Evi_j(j) = J_ini_; % History of resulting neg. log-evidence
                    Xout_ini = Xout_ini_;
                else
                    % Storing values to array
                     Xout_j{j} = Xout_i(:,:,i_min);
                        w_j{j} = w_i(:,:,i_min); % History of parameters of terms selected
                    J_Evi_j(j) = J_Evi_i(i_min); % History of resulting neg. log-evidence
                    Xout_ini = Xout_i(:,:,i_min);
                end
    
                % Visualise
                if options.PlotXout
                    visualise_denoise(Xout_ini);
                end
                if options.VerboseLevel >=1
                    disp(['After removal of the ' num2str(j) 'th term']);
                    disp(w_j{j});
                end
    
                if j > 2
                    if J_Evi_j(j)>=J_Evi_j(j-1)
                        ContRun_flag=ContRun_flag+1; % Raise the counter if J_evi start rising
                    else
                        ContRun_flag=0; % Reset counter if J_evi drops. i.e. Only if J_Evi_j rise continuously will the run exit
                    end
                end
    
                mask=mask_;
                j=j+1;
            else
                % If decided not to proceed, retry the run with another
                % Xout_ini
                Xout_ini = Xout_ini_; 
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
end

function [proceed_flag,mask_,i_min]=RemovalDecision(mask,J_Evi_i,m,ND)
    if all(J_Evi_i==Inf)
        mask_= mask; % Run has failed. Re-run the previous initial trial to get better Xout_ini
        proceed_flag = false;
        i_min=0;
    else
        % Select the index with the MINIMUM J (maximum posterior likelihood).
        % This means: "when you drop the polynomial term with this index, the remaining model
        % has higher evidence than all the other alternative models."
        [~,i_sorted]=sort(J_Evi_i);
    
        % Remove this index from the indices that can be selected next
        % (avoid making no active term in each column) 
    
        iii=0;
        M= numel(mask)/ND;
        mask_next = false(M,ND);
        while ~all(any(mask_next,1)) % Keep going if any column has no true term
            iii=iii+1;
            mask_next = reshape(mask,[],ND);
            mask_next(m(i_sorted(iii))) = false;           
        end
        mask_ = reshape(mask_next,M*ND,1);
        i_min = i_sorted(iii);
    
        % Reset Initial Trial Call
        if J_Evi_i(i_min)==Inf % In case if all Inf except the only term that is also single on the
        % column
            proceed_flag=false;
            mask_= mask; 
            i_min=0;
        else
            proceed_flag = true;
        end
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