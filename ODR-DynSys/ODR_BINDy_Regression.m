function [w,J_Evi,Xout,exitflag]=ODR_BINDy_Regression(data,libs,mask_loc,lsqnonlinopt,InitGuessOpt,options,varargin)
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
        if options.LinUseDenoise 
            plin = EnsembleLinRegress_Wprior(reshape(x_ini,data.Nx,data.ND),...
                data.Theta_fun,data.DMat,data.IMat,mask_loc,data.SigmaY,data.SigmaBeta,InitGuessOpt);
        else
            plin = EnsembleLinRegress_Wprior(data.Xdata,...
                data.Theta_fun,data.DMat,data.IMat,mask_loc,data.SigmaY,data.SigmaBeta,InitGuessOpt);
        end
        
        try
        targetfun = @(Xp)Jsq_xp(Xp,data,mask_loc);
        [Xp,resnorm,resf,exitflag,output] = lsqnonlin(targetfun,...
            [x_ini;reshape(plin,[],1)],...
            [],[],lsqnonlinopt);
        
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

        catch
            if options.VerboseLevel >=2
                warning('Something is wrong in lsqnonlin');
            end
            exitflag=0;
            w = zeros(libs.M,data.ND);
            Xout = zeros(data.Nx,data.ND);
            J_Evi = Inf;
        end

        if options.VerboseLevel >=2
            if exitflag<=0
                warning(['lsqnonlin cannot converge:' num2str(exitflag)]);
            elseif any(eig(Hessian_p)<0)
                warning(['Hessian not positive-def. Exitflag:' num2str(exitflag) ]);
            end
        end
        
        jjj = jjj +1;
    end
    if exitflag<=0
        % Do not remove term if convergence failed.
        J_Evi= Inf;
    elseif jjj >= options.MaxNoiseIter
        if options.VerboseLevel >=2
            warning('Over Iterations.');
        end
        J_Evi= Inf;
        % J_Evi= resnorm/2 ...
        %     +(log(2*pi)-log(alpha))*Np/2; ... ; % WIP: NOT SURE if it's right
    % elseif exitflag <=0
    %     warning('lsqnonlin did not converge');
    %     J_Evi = Inf;
    else
        try 
            J_logdet_Hessian = sum(log(diag(chol(Hessian_p/2/pi))));
        catch
            J_logdet_Hessian = logdet(Hessian_p/2/pi)/2;
        end
        % Negative Log-Evidence, by Laplace Approximation WIP!!
        J_Evi=resnorm/2 ...
            ...%+log(2*pi)*Neq/2+sum(log(CoV(2,2))/2)... % Constant
            ...%+log(2*pi)*(ND*Neq)/2+sum(log(CoV(1,1))/2,'all')... % Constant
            +(log(2*pi))*sum(mask_loc,'all')/2 +sum(log(data.SigmaBeta(mask_loc)))... % Change with # param % TODO: REMOVING THIS LINE IMPROVE THINGS??
            +J_logdet_Hessian; % Need to check sign
    end
end