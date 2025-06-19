function out = Polynomial4D3O(X)
N = numel(X)/4;
X = reshape(X,N,4);
   x1=X(:,1);
   x2=X(:,2);
   x3=X(:,3);
   x4=X(:,4);
out = [ones(N,1) x1 x2 x3 x4 x1.^2 x1.*x2 x1.*x3 x1.*x4 x2.^2 x2.*x3 x2.*x4 x3.^2 x3.*x4 x4.^2 x1.^3 x1.^2.*x2 x1.^2.*x3 x1.^2.*x4 x1.*x2.^2 x1.*x2.*x3 x1.*x2.*x4 x1.*x3.^2 x1.*x3.*x4 x1.*x4.^2 x2.^3 x2.^2.*x3 x2.^2.*x4 x2.*x3.^2 x2.*x3.*x4 x2.*x4.^2 x3.^3 x3.^2.*x4 x3.*x4.^2 x4.^3 ];
end