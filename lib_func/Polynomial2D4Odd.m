function out = Polynomial2D4Odd(X)
N = numel(X)/2;
X = reshape(X,N,2);
   x1=X(:,1);
   x2=X(:,2);
out = zeros(N,15,2,2);
out(:,4,1,1) = 2;
out(:,7,1,1) = 6.*x1;
out(:,8,1,1) = 2.*x2;
out(:,11,1,1) = 12.*x1.^2;
out(:,12,1,1) = 6.*x1.*x2;
out(:,13,1,1) = 2.*x2.^2;
out(:,5,1,2) = 1;
out(:,8,1,2) = 2.*x1;
out(:,9,1,2) = 2.*x2;
out(:,12,1,2) = 3.*x1.^2;
out(:,13,1,2) = 4.*x1.*x2;
out(:,14,1,2) = 3.*x2.^2;
out(:,5,2,1) = 1;
out(:,8,2,1) = 2.*x1;
out(:,9,2,1) = 2.*x2;
out(:,12,2,1) = 3.*x1.^2;
out(:,13,2,1) = 4.*x1.*x2;
out(:,14,2,1) = 3.*x2.^2;
out(:,6,2,2) = 2;
out(:,9,2,2) = 2.*x1;
out(:,10,2,2) = 6.*x2;
out(:,13,2,2) = 2.*x1.^2;
out(:,14,2,2) = 6.*x1.*x2;
out(:,15,2,2) = 12.*x2.^2;
end