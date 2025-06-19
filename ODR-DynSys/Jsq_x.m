function [out,J] = Jsq_x(x,p,Data,mask)
    %% Elements that make up a Gaussian Likelihood J(x) given parameter p
    %  assuming a 2-norm (Gaussian) log-loss J = 1/2 (delta^2 + eta^2),
    %  where delta = xdot-Theta(x)*p / SigmaY , eta = x-xhat / SigmaX
    %  xdot estimated by DMat*x and Theta(x) = IMAT*lib(x)*p
    %  xhat being the data of x
    %  This function gives the breakdown of the J in its elements, 
    %  i.e. [eta;delta]
    %  The 2-norm of this vector gives J*2.
    %  
    %  OUTPUT:
    %     out : [eta;delta]
    %
    %  INPUT:
    %     x: Local guess value of x we are regressing
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

    delta = zeros(Data.Nx,Data.ND);
    eta = zeros(Data.Neq,Data.ND);
    if nargout >1
        dD = Data.dTheta_fun(x);
        df = NaN(Data.Nx,Data.ND,Data.ND);
    end
    for d = 1:Data.ND
        startd=sum(mask(:,1:d-1),'all')+1;
        endd = startd-1+sum(mask(:,d));
        plocd = p(startd:endd);
        delta(:,d) = (x(:,d)-Data.Xdata(:,d))./Data.SigmaX(:,d);
        eta(:,d) = (dx(:,d)-Theta_loc(:,mask(:,d))*plocd)./Data.SigmaY(:,d);
        if nargout >1
            df(:,d,:) = pagemtimes(dD(:,mask(:,d),:),plocd);
        end
    end
    out = [reshape(eta,Data.Neq*Data.ND,1);reshape(delta,Data.Nx*Data.ND,1)];
    if nargout >1
        detadx = (kron(speye(Data.ND),Data.DMat)...
        -kron(speye(Data.ND),Data.IMat)*...
        sparse(repmat((1:Data.Nx*Data.ND)',Data.ND,1),...
        reshape(repmat((1:Data.Nx)',Data.ND,Data.ND)+(0:Data.ND-1)*Data.Nx,Data.Nx*Data.ND*Data.ND,1),...
        reshape(df,[],1),Data.Nx*Data.ND,Data.Nx*Data.ND))...
        ./reshape(Data.SigmaY,[],1);
        J = [detadx;
            spdiags(1./Data.SigmaX(:),0,Data.Nx*Data.ND,Data.Nx*Data.ND)];
    end
end