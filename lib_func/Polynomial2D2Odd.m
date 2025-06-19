function out = Polynomial2D2Odd(X)
N = numel(X)/2;
X = reshape(X,N,2);
   x1=X(:,1);
   x2=X(:,2);
out = zeros(N,6,2,2);
out(:,4,1,1) = 2;
out(:,5,1,2) = 1;
out(:,5,2,1) = 1;
out(:,6,2,2) = 2;
end