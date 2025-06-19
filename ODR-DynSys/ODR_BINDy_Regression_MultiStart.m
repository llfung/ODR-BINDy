function [w,J_Evi,Xout,exitflag]=ODR_BINDy_Regression_MultiStart(data,libs,mask_loc,lsqnonlinopt,InitGuessOpt,options,varargin)
    % Copyright 2025, All Rights Reserved
    % Code by Lloyd Fung
    %    For package ODR-BINDy
    if nargin > 6
        x_ini = reshape(varargin{1},data.Nx*data.ND,1);
    else
        x_ini = reshape(data.Xdata,data.Nx*data.ND,1);
    end
    jjj=0;
    Hessian_p = -1;
    exitflag = 1;
    while (any(eig(Hessian_p)<0) && exitflag>0) && jjj < options.MaxNoiseIter
        % To get an initial guess for p, we perform ensemble linear regression  
        plin_array = NaN(options.StartPtNum,sum(mask_loc,"all"));
        if options.LinUseDenoise
            parfor kk = 1:options.StartPtNum
            plin_array(kk,:) = (EnsembleLinRegress_Wprior(reshape(x_ini,data.Nx,data.ND),...
                data.Theta_fun,data.DMat,data.IMat,mask_loc,data.SigmaY,data.SigmaBeta,InitGuessOpt))';
            end
        else
            parfor kk = 1:options.StartPtNum
            plin_array(kk,:) = (EnsembleLinRegress_Wprior(data.Xdata,...
                data.Theta_fun,data.DMat,data.IMat,mask_loc,data.SigmaY,data.SigmaBeta,InitGuessOpt))';
            end
        end
        tpoints = CustomStartPointSet(...
            [ones(options.StartPtNum,1)*x_ini' plin_array]);

        ms = MultiStart('Display','iter','UseParallel',true,...
            'FunctionTolerance',1e-8,'XTolerance',options.StepTolerance);

        if options.VerboseLevel >=3
            ms.Display = 'iter';
        elseif options.VerboseLevel >=2
            ms.Display = 'final';
        else
            ms.Display = 'off';
        end

        targetfun = @(Xp)Jsq_xp(Xp,data,mask_loc);
        problem = createOptimProblem('lsqnonlin','objective',targetfun,...
            'options',lsqnonlinopt,'x0',[x_ini;plin_array(1,:)']);


        [Xp,resnorm,exitflag,output] = run(ms,problem,tpoints);
        
        % disp(['Iterations:' num2str(output.iterations)]);
        xtemp = Xp(1:data.Nx*data.ND);
        ptemp = Xp(data.Nx*data.ND+1:end,1);
        % [xtemp,xnorm] = lsqnonlin(@(x)numcheck(x,ptemp,data,mask_loc),...
        %     xtemp,[],[],lsqnonlinopt); % Final correction to improve accruacy
        
        
        if exitflag>0
            if options.UseGaussNewtonEst
                Hessian_p = d2Jdp2_GaussNewtonEst(reshape(xtemp,data.Nx,data.ND),ptemp,libs,data,mask_loc);
            else
            % Hessian_p = d2Jdp2(reshape(xtemp,data.Nx,data.ND),ptemp,libs,data,mask_loc);
            % if any(eig(Hessian_p)<0)
                
                xtemp = lsqnonlin(@(x)Jsq_x(x,ptemp,data,mask_loc),xtemp,[],[],lsqnonlinopt);
                xtemp = lsqnonlin(@(x)dJdx(x,ptemp,data,mask_loc,libs),xtemp,[],[],lsqnonlinopt);
                
                Hessian_p = d2Jdp2(reshape(xtemp,data.Nx,data.ND),ptemp,libs,data,mask_loc);
            % end
            end
        end
        Xout = reshape(xtemp,data.Nx,data.ND);
        w = zeros(libs.M,data.ND);
        w(mask_loc) = Xp(data.Nx*data.ND+1:end,1);
        
        % Calculate Hessian
        % Hessian_p = d2Jdp2(Xout,ptemp,libs,data,mask_loc);

        % Visualise
        % figure(1);
        % plot(data.Xdata,'x');hold on;
        % plot(Xout);hold off;
        
        if exitflag<=0
            warning(['lsqnonlin cannot converge:' num2str(exitflag)]);
        elseif any(eig(Hessian_p)<0)
            warning(['Hessian not positive-def. Exitflag:' num2str(exitflag) ]);
        end
        
        jjj = jjj +1;
    end
    if exitflag<=0
        % Do not remove term if convergence failed.
        J_Evi= Inf;
    elseif jjj >= options.MaxNoiseIter
        warning('Over Iterations.');
        J_Evi= Inf;
        % J_Evi= resnorm/2 ...
        %     +(log(2*pi)-log(alpha))*Np/2; ... ; % WIP: NOT SURE if it's right
    % elseif exitflag <=0
    %     warning('lsqnonlin did not converge');
    %     J_Evi = Inf;
    else
        % Negative Log-Evidence, by Laplace Approximation WIP!!
        J_Evi=resnorm/2 ...
            ...%+log(2*pi)*Neq/2+sum(log(CoV(2,2))/2)... % Constant
            ...%+log(2*pi)*(ND*Neq)/2+sum(log(CoV(1,1))/2,'all')... % Constant
            +(log(2*pi))*sum(mask_loc,'all')/2 +sum(log(data.SigmaBeta(mask_loc)))... % Change with # param % TODO: REMOVING THIS LINE IMPROVE THINGS??
            +logdet(Hessian_p/2/pi)/2; % Need to check sign
    end
end