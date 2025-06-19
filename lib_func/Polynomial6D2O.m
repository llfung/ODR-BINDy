function out = Polynomial6D2O(X)
N = numel(X)/6;
X = reshape(X,N,6);
   x1=X(:,1);
   x2=X(:,2);
   x3=X(:,3);
   x4=X(:,4);
   x5=X(:,5);
   x6=X(:,6);
out = [ones(N,1) x1 x2 x3 x4 x5 x6 x1.^2 x1.*x2 x1.*x3 x1.*x4 x1.*x5 x1.*x6 x2.^2 x2.*x3 x2.*x4 x2.*x5 x2.*x6 x3.^2 x3.*x4 x3.*x5 x3.*x6 x4.^2 x4.*x5 x4.*x6 x5.^2 x5.*x6 x6.^2 ];
end