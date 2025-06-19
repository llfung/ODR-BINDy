function [out,Jacobian] = dJdx(x,p,Data,mask,libs)
    %% Gradient of Gaussian loss J against x given parameter p
    %  assuming a 2-norm (Gaussian) log-loss J = 1/2 (delta^2/2 + eta^2),
    %  where delta = xdot-Theta(x)*p / SigmaY , eta = x-xhat / SigmaX
    %  xdot estimated by DMat*x and Theta(x) = IMAT*lib(x)*p
    %  xhat being the data of x
    %  
    %  OUTPUT:
    %     out : dJdx at give p
    %     Jacobian: Jacobian of dJdx against x (i.e. Hessian!). For
    %     numerical convergence
    %
    %  INPUT:
    %     X: Local guess value of X we are regressing
    %     p: Given value of parameter
    %     Data: Struct containing the necessary data
    %        XData: The observed X (Nx x ND)
    %        DMat : (Sparse) matrix to get estimate of dx (Neq x Nx)
    %        IMat : (Sparse) matrix to map local value of Theta(x)*p 
    %               to corresponds to dx (usually identity if FD is used, Neq x Nx)
    %        SigmaY: Noise STD in the dx dimension (Neq x ND)
    %        SigmaX: Noise STD in the x dimension (Nx x ND)
    %     mask: Parameter library mask (sum(mask,'all') = length(p), M x ND)
    %
    % Copyright 2025, All Rights Reserved
    % Code by Lloyd Fung
    %    For package ODR-BINDy
    %


    mask = reshape(mask,[],Data.ND);
    x = reshape(x(1:Data.Nx*Data.ND),Data.Nx,Data.ND);
    Theta_loc = Data.IMat*Data.Theta_fun(x);
    dx = Data.DMat*x;
    dD = libs.dTheta_fun(x);
    if nargout > 1
        ddD = libs.ddTheta_fun(x);
    end

    delta = zeros(Data.Nx,Data.ND);
    eta = zeros(Data.Neq,Data.ND);
    df = zeros(Data.Nx,Data.ND,Data.ND);
    if nargout > 1
        ddf = zeros(Data.Nx,Data.ND,Data.ND,Data.ND);
    end
    for d = 1:Data.ND
        startd=sum(mask(:,1:d-1),'all')+1;
        endd = startd-1+sum(mask(:,d));
        p_d = p(startd:endd);
        delta(:,d) = (x(:,d)-Data.Xdata(:,d))./Data.SigmaX(:,d).^2;
        eta(:,d) = (dx(:,d)-Theta_loc(:,mask(:,d))*p_d)./Data.SigmaY(:,d).^2;
        df(:,d,:) = pagemtimes(dD(:,mask(:,d),:),p_d);
        if nargout > 1
            ddf(:,d,:,:) = pagemtimes(ddD(:,mask(:,d),:,:),p_d);
        end
    end
    df_diag = sparse(repmat((1:Data.Nx*Data.ND)',Data.ND,1),...
        reshape(repmat((1:Data.Nx)',Data.ND,Data.ND)+(0:Data.ND-1)*Data.Nx,Data.Nx*Data.ND*Data.ND,1),...
        reshape(df,[],1),Data.Nx*Data.ND,Data.Nx*Data.ND);
    deta = (kron(speye(Data.ND),Data.DMat)...
        -kron(speye(Data.ND),Data.IMat)*df_diag);

    out = (reshape(eta,1,Data.Neq*Data.ND)*deta)'+reshape(delta,Data.Nx*Data.ND,1);

    if nargout > 1
       deta_ = deta./reshape(Data.SigmaY,[],1);
    
        % Find first derivative dx/dbeta
        eta_i_eta_iab= permute(...
            pagemtimes(permute((-Data.IMat')*(eta./Data.SigmaY.^2),[3,2,1]),...
            permute(ddf,[2,5,1,3,4])),...
            [3,4,5,1,2]);
        eta_i_eta_iab_diag = sparse(repmat((1:Data.Nx*Data.ND)',Data.ND,1),...
            reshape(repmat((1:Data.Nx)',Data.ND,Data.ND)+(0:Data.ND-1)*Data.Nx,Data.Nx*Data.ND*Data.ND,1),...
            reshape(eta_i_eta_iab,[],1),Data.Nx*Data.ND,Data.Nx*Data.ND);
    
        Jacobian = deta_'*deta_ + eta_i_eta_iab_diag + spdiags(reshape(1./Data.SigmaX.^2,[],1),0,Data.Nx*Data.ND,Data.Nx*Data.ND);
    end
end
