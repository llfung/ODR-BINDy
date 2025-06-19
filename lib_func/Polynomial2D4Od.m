function out = Polynomial2D4Od(X)
N = numel(X)/2;
X = reshape(X,N,2);
   x1=X(:,1);
   x2=X(:,2);
out = zeros(N,15,2);
out(:,:,1) = [zeros(N,1) ones(N,1) zeros(N,1) 2.*x1 x2 zeros(N,1) 3.*x1.^2 2.*x1.*x2 x2.^2 zeros(N,1) 4.*x1.^3 3.*x1.^2.*x2 2.*x1.*x2.^2 x2.^3 zeros(N,1) ];
out(:,:,2) = [zeros(N,1) zeros(N,1) ones(N,1) zeros(N,1) x1 2.*x2 zeros(N,1) x1.^2 2.*x1.*x2 3.*x2.^2 zeros(N,1) x1.^3 2.*x1.^2.*x2 3.*x1.*x2.^2 4.*x2.^3 ];
end