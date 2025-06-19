function out = d2Jdp2_GaussNewtonEst(x,p,libs,Data,mask)
    %% Hessian of log-posterior J(x,p) against p given x has converged to minimum given p
    %  assuming a Gaussian log-posterior J = 1/2 (delta^2/2 + eta^2 + dp^2),
    %  where Likelihood: delta = xdot-Theta(x)*p / SigmaY , eta = x-xhat / SigmaX
    %    and Prior: dp = p / SigmaBeta (0 mean prior with SigmaBeta^2 variance)
    %  xdot estimated by DMat*x and Theta(x) = IMAT*lib(x)*p
    %  xhat being the data of x
    %  This function gives the Hessian of Jp(p)=J(x_min,p), 
    %  given x_min has converged to the minimum of J(x,p) at given p.
    %  
    %  OUTPUT:
    %     out : dJdx at give p
    %     Jacobian: Jacobian of dJdx against x (i.e. Hessian!). For
    %     numerical convergence
    %
    %  INPUT:
    %     x: Value of x_min (ASSUMED fulfill dJdx(x_min) = 0 !!)
    %     p: Local value of p 
    %     Data: Struct containing the necessary data
    %        XData: The observed X (Nx x ND)
    %        DMat : (Sparse) matrix to get estimate of dx (Neq x Nx)
    %        IMat : (Sparse) matrix to map local value of Theta(x)*p 
    %               to corresponds to dx (usually identity if FD is used, Neq x Nx)
    %        SigmaY: Noise STD in the dx dimension (Neq x ND)
    %        SigmaX: Noise STD in the x dimension (Nx x ND)
    %     mask: Parameter library mask (sum(mask,'all') = length(p), M x ND)
    %     libs: Library of Candidate functions and their derivatives w.r.t. x
    %
    %  N.B.
    %     If p=p_min that minimise J also, then eig(d2Jdp2) should have all
    %     positive eigenvalues
    %     If not, it's likely some convergence issue with x or p.
    %
    % Copyright 2025, All Rights Reserved
    % Code by Lloyd Fung
    %    For package ODR-BINDy
    %

    mask = reshape(mask,[],Data.ND);
    Nx = Data.Nx;
    Neq = Data.Neq;
    ND = Data.ND;
    Np = length(p);

    D = libs.Theta_fun(x);
    dD = libs.dTheta_fun(x);
    ddD = libs.ddTheta_fun(x);

    f = NaN(Nx,ND);
    df = NaN(Nx,ND,ND);
    ddf = NaN(Nx,ND,ND,ND);
    % dddf = NaN(Nx,ND,ND,ND,ND);

    pd = zeros(Np,1);
    pd_ind = cell(ND,1);
    for d=1:ND
        startd=sum(mask(:,1:d-1),'all')+1;
        endd = startd-1+sum(mask(:,d));
        ps = startd:endd;
        pd(ps) = d;
        pd_ind{d} = ps;
        plocd = p(ps);

        f(:,d) = D(:,mask(:,d))*plocd;
        df(:,d,:) = pagemtimes(dD(:,mask(:,d),:),plocd);
        ddf(:,d,:,:) = pagemtimes(ddD(:,mask(:,d),:,:),plocd);
        % dddf(:,d,:,:,:) = reshape(libs.dddTheta_fun_f(x,plocd,mask(:,d)),Nx,1,ND,ND,ND);
    end

    eta = Data.DMat*reshape(x,[],ND)-Data.IMat*f;
    
    %% %%%%%%%%%%%%%%% dxdp  &  d2xdp2  %%%%%%%%%%%%%%%%%%%%%
    % eta_i,b (with IMat) (Neq*ND,Nx*ND)
    df_diag = sparse(repmat((1:Nx*ND)',ND,1),...
        reshape(repmat((1:Nx)',ND,ND)+(0:ND-1)*Nx,Nx*ND*ND,1),...
        reshape(df,[],1),Nx*ND,Nx*ND);
    eta_ib = (kron(speye(ND),Data.DMat)...
        -kron(speye(ND),Data.IMat)*df_diag)./reshape(Data.SigmaY,[],1);

    % eta_i,m
    eta_im = zeros(Nx,ND,Np);
    for d = 1:ND
        eta_im(:,d,pd_ind{d}) = D(:,mask(:,d));
    end
    
    % eta_i,bm (i,d,m,b)
    eta_imb = zeros(Nx,ND,Np,ND);
    for d = 1:ND
        eta_imb(:,d,pd_ind{d},:) = reshape(dD(:,mask(:,d),:),Nx,1,length(pd_ind{d}),ND);
    end

    % Find first derivative dx/dp
    eta_i_eta_iab= permute(...
        pagemtimes(permute((-Data.IMat')*(eta./Data.SigmaY.^2),[3,2,1]),...
        permute(ddf,[2,5,1,3,4])),...
        [3,4,5,1,2]);
    eta_i_eta_iab_diag = sparse(repmat((1:Nx*ND)',ND,1),...
        reshape(repmat((1:Nx)',ND,ND)+(0:ND-1)*Nx,Nx*ND*ND,1),...
        reshape(eta_i_eta_iab,[],1),Nx*ND,Nx*ND);

    LHS = eta_ib'*eta_ib + eta_i_eta_iab_diag + spdiags(reshape(1./Data.SigmaX.^2,[],1),0,Nx*ND,Nx*ND);

    eta_i_eta_ian = permute(...
        pagemtimes(permute((-Data.IMat)'*(eta./Data.SigmaY.^2),[3,2,1]),...
        permute(eta_imb,[2,3,1,4])),...
        [3,4,2,1]);
    RHS = -reshape(eta_i_eta_ian,Nx*ND,Np) ...
        -eta_ib'*reshape(reshape((-Data.IMat)*reshape(eta_im,Nx,ND*Np),Neq,ND,Np)./Data.SigmaY,Neq*ND,Np);

    dxdp = LHS \ RHS;

    % IMat and sigmaY is only added during the final i contraction 
    dxdp_ = permute(reshape(dxdp,Nx,ND,Np),[2,3,1]); %(i x b,n) -> (b,n,i)

    % eta_i,ab * dxdbeat (a,n,i,d) -> (i,d,a,n)
    % eta_iab_dxbdbeta  = permute(pagemtimes(permute(ddf,[3,4,1,2]),...
    %     dxdp_),[3,4,1,2]);

    % Start finding second derivatives d2x/dbeta2
    % c11 = permute(pagemtimes(...
    %     permute(reshape(-Data.IMat'*reshape((reshape(Data.SigmaY.^-1,[],1).*eta_ib)*dxdp,Neq,ND*Np),Nx,ND,Np),[3,2,1]),... %eta_i,b * dxdp (n,d,i)
    %     permute(eta_imb,[2,3,1,4])),... %eta_i,ma (d,m,i,a)
    %         [3,4,1,2]);
    % 
    % c12 = reshape(...
    %     (eta_ib')*(reshape(Data.SigmaY.^-1,[],1).*kron(speye(ND),-Data.IMat)*...
    %     reshape(permute(...
    %     pagemtimes(permute(eta_imb,[2,4,1,3]),dxdp_),... % eta_imb * dxdp (d,n,i,m)
    %     [3,1,2,4]),Nx*ND,Np*Np)),...
    %     Nx*ND,Np,Np); % eta_imb * dxdp (i x d,n,m)
    % 
    % c13 = permute(pagemtimes( ...
    %     permute(...
    %     reshape((-Data.IMat)'*reshape(...
    %     reshape(-Data.IMat*reshape(eta_im,Nx,ND*Np),Neq,ND,Np)./(Data.SigmaY.^2),...
    %     Neq,ND*Np),Nx,ND,Np),... % eta_i,m (j,d,m)
    %     [3,2,1]),... % eta_i,m (j,d,m)-> (m,d,j)
    %     permute(eta_iab_dxbdbeta,[2,3,1,4])),... %eta_i,ab * dxdp (j,d,a,n)->(d,a,j,n)
    %     [3,2,4,1]); %(m,a,j,n)-> (j,a,n,m)
    % 
    % c14 = zeros(Nx,ND,Np,Np);
    % for d0 = 1:ND
    %     c14(:,:,pd_ind{d0},:) = ((-Data.IMat)'*(eta(:,d0)./Data.SigmaY(:,d0).^2)).*...
    %                 permute(...
    %                 pagemtimes(permute(ddD(:,mask(:,d0),:,:),[2,3,1,4]),... %ddD (i,m,a,b) -> (m,a,i,b)
    %                 dxdp_),... %dxdp (i,a,n) -> (a,n,i)
    %                 [3,4,1,2]); %(m,n,i,b)->(i,b,m,n)
    % end
    % 
    % c31 = reshape(permute(pagemtimes( ...
    %     permute(reshape(-Data.IMat'*reshape((reshape(Data.SigmaY.^-1,[],1).*eta_ib)*dxdp,Neq,ND*Np),Nx,ND,Np),[3,2,1]),... % eta_i,c dxdp (n,d,i)
    %     permute(eta_iab_dxbdbeta,[2,4,1,3])),... %eta_i,ab * dxdp (d,m,i,a)
    %     [3,4,1,2]),Nx*ND,Np,Np);
    % 
    % c32 = eta_ib'*(reshape(Data.SigmaY.^-1,[],1).*kron(speye(ND),-Data.IMat)*...
    %     reshape(permute(pagemtimes(...
    %     permute(dxdp_,[2,1,3]),... %dxdp (n,a,i)
    %     pagemtimes(...
    %     permute(ddf,[3,4,1,2]),...% ddf (i,d,a,b) - > (a,b,i,d)
    %     dxdp_)), ... %dxdp (b,m,i)
    %     [3,4,1,2]),Nx*ND,Np*Np)); %  dxdp * ddf * dxdp (n,m,i,d)-> (i,d,n,m) -> (i x d, n x m)
    % 
    % c34 = permute(pagemtimes(permute((-Data.IMat)'*(eta./Data.SigmaY.^2),[3,2,1]),... %(i,d) -> (_,d,i)
    %         permute(...
    %         pagemtimes(permute(dxdp_,[2,1,3]),... %dxdp (i,c,m) -> (c,m,i) -> (m,c,i)
    %         pagemtimes(permute(dddf,[3,4,1,2,5]),... %ddD (i,d,c,b,a) -> (c,b,i,d,a)
    %         dxdp_)),... %dxdp (i,b,n) -> (b,n,i)
    %         [4,5,3,1,2])),...  %(m,n,i,d,a)->(d,a,i,m,n)
    %         [3,2,4,5,1]);
    % 
    % c42 = reshape(...
    %         permute(...
    %             pagemtimes(...
    %                 permute(reshape(...
    %                     Data.IMat'*reshape(...
    %                         reshape(Data.IMat*reshape(eta_im,Nx,ND*Np),Neq,ND,Np)./Data.SigmaY.^2,...
    %                     Neq,ND*Np),...
    %                 Nx,ND,Np),[3,2,1]),... %eta_i,m (i, d, n) ->(n,d,i)
    %             permute(eta_imb,[2,3,1,4])),... % eta_imb (i,d,m,b) -> (d,m,i,b)
    %         [3,4,1,2]),...
    %       Nx*ND,Np,Np); %(n,m,i,b)->(i,b,n,m)
    % 
    % % a,n,m
    % c1temp = reshape(c11(:) + c12(:) + c13(:) + c14(:) ,Nx*ND,Np,Np);
    % c3temp = c31 + permute(c31,[1,3,2]);
    % c3temp = c3temp + reshape(c32,Nx*ND,Np,Np) + reshape(c34,Nx*ND,Np,Np);
    % c4temp = c42 + permute(c42,[1,3,2]);
    % 
    % d2xdp2= reshape(...
    %     LHS \ reshape(-c1temp - permute(c1temp,[1,3,2]) - c3temp - c4temp,Nx*ND,Np*Np),...
    %     Nx*ND,Np,Np); % check and consistent with slow version

    %% %%%%% Calculating d2Jdp2  %%%%%%%%
    % eta_temp = reshape(reshape(eta./Data.SigmaY.^2,1,[])*kron(speye(ND),-Data.IMat)*...
    %     reshape(permute(...
    %     pagemtimes(permute(eta_imb,[2,4,1,3]),dxdp_),... % eta_imb * dxdp (d,n,i,m)
    %     [3,1,2,4]),Nx*ND,Np*Np),Np,Np);
    % eta_temp = eta_temp + eta_temp';
    % eta_temp = eta_temp + reshape(reshape(eta./Data.SigmaY,1,[])*...
    %     eta_ib*reshape(d2xdp2,Nx*ND,Np*Np),Np,Np);
    eta_temp2 = eta_ib*dxdp+...
        reshape(reshape(-Data.IMat*reshape(eta_im,Nx,ND*Np),Neq,ND,Np)./Data.SigmaY,Neq*ND,Np);
    % d2J_eta = eta_temp + eta_temp2'*eta_temp2; 
    d2J_eta = eta_temp2'*eta_temp2;
    % d2J_x = tensorprod(d2xdp2,reshape((x-Data.Xdata)./Data.SigmaX.^2,[],1),1,1)...
    d2J_x =  dxdp'*diag(reshape(Data.SigmaX.^-2,[],1))*dxdp;
    out = d2J_eta+d2J_x+diag(Data.SigmaBeta(mask).^-2);

    %% %%%%%%%%%%%%%%%%%%%%%%%%%% Numerical CHECKS %%%%%%%%%%%%%%%%%%%
    % [c11_check,c12_check,c13_check,c31_check,c32_check,c42_check] = eta2_numcheck(Data,mask,x,p,dxdp);
    % [d2xdp11,d2xdp12,d2xdp22] = d2xdp2_numcheck(Data,mask,x,p,6,7);

    % if det(d2Jdp2)<0
    % d2Jdp2_ = NaN(Np,Np);
    % for p_select1 = 1:Np
    %     for p_select2 = 1:Np
    %         if p_select1~=p_select2
    %             [d2p11,d2p12] = dx_numcheck(Data,mask,x,p,p_select1,p_select2);
    %             d2Jdp2_(p_select1,p_select1) = d2p11;
    %             d2Jdp2_(p_select1,p_select2) = d2p12;
    %         end
    %     end
    % end
    % disp(det(d2Jdp2_));
    % end
    % %%%%%%%%%%%%%%%%%%%%%%%%%%% Numerical CHECKS %%%%%%%%%%%%%%%%%%%

end

