function wlin = EnsembleLinRegress_Wprior(Xdata_lin,Theta_fun,DMat,IMat,mask_loc,SigmaY,SigmaBeta,options)
    %% Initial Parameter guess by Ensembling Regression
    % Used in ODR-BINDy as a parameter p initial condition guesser
    %
    % Copyright 2025, All Rights Reserved
    % Code by Lloyd Fung
    %    For package ODR-BINDy

    ND = size(Xdata_lin,2);
    Neq = size(IMat,1);
    M = numel(mask_loc)/ND;
    wlin = [];

    Theta_all = IMat*Theta_fun(Xdata_lin);
    
    for d = 1:ND
        % Select the Theta column according to indices selected (mark_loc)
        mask_locd = mask_loc((d-1)*M+1:d*M,1);
        SigmaBeta_loc = SigmaBeta(mask_locd);
        
        Theta = Theta_all(:,mask_locd);
        dx = (DMat*Xdata_lin(:,d))./SigmaY(:,d).^2; % TODO: Need better est. of what the noise variance should be. How about SigmaX?

        % Find parameter value by ensembled linear regression! (Prior serve as regularisation)
        if options.DataSelectRatio==1
            % Proper Bootstrapping
            coeff= bootstrp(options.EnsembleSampleSize,@(Theta_,dx,SigmaY_)(Theta_'*(Theta_./SigmaY_.^2)+diag(SigmaBeta_loc.^-2))... % LHS
                \(Theta_'*dx),Theta,dx,SigmaY(:,d));
            if options.Bragging
                coeff_bragging = median(coeff,1)';
                wlin = [wlin;coeff_bragging];
            else
                coeff_bagging = mean(coeff,1)';
                wlin = [wlin;coeff_bagging];
            end
            
        else
            % Bootstrapping by sampling slightly less
            coeff_array = zeros(size(Theta,2),options.EnsembleSampleSize);
            for ll=1:options.EnsembleSampleSize
                ind = randi(Neq,floor(Neq*DataSelectRatio),1); % Sampling with replacement
                Theta_ = Theta(ind,:);
                coeff_array(:,ll) =  (Theta_'*(Theta_./SigmaY(ind,d).^2)+diag(SigmaBeta_loc.^-2))... % LHS
                    \(Theta_'*dx(ind,:)); % RHS
            end
            if options.Bragging
                wlin = [wlin;median(coeff_array,2)];
            else
                wlin = [wlin;mean(coeff_array,2)]; 
            end
        end
    end
end