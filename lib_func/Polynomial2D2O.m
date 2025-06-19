function out = Polynomial2D2O(X)
N = numel(X)/2;
X = reshape(X,N,2);
   x1=X(:,1);
   x2=X(:,2);
out = [ones(N,1) x1 x2 x1.^2 x1.*x2 x2.^2 ];
end