function out = Polynomial5D2Odd(X)
N = numel(X)/5;
% X = reshape(X,N,5);
   % x1=X(:,1);
   % x2=X(:,2);
   % x3=X(:,3);
   % x4=X(:,4);
   % x5=X(:,5);
out = zeros(N,21,5,5);
out(:,7,1,1) = 2;
out(:,8,1,2) = 1;
out(:,9,1,3) = 1;
out(:,10,1,4) = 1;
out(:,11,1,5) = 1;
out(:,8,2,1) = 1;
out(:,12,2,2) = 2;
out(:,13,2,3) = 1;
out(:,14,2,4) = 1;
out(:,15,2,5) = 1;
out(:,9,3,1) = 1;
out(:,13,3,2) = 1;
out(:,16,3,3) = 2;
out(:,17,3,4) = 1;
out(:,18,3,5) = 1;
out(:,10,4,1) = 1;
out(:,14,4,2) = 1;
out(:,17,4,3) = 1;
out(:,19,4,4) = 2;
out(:,20,4,5) = 1;
out(:,11,5,1) = 1;
out(:,15,5,2) = 1;
out(:,18,5,3) = 1;
out(:,20,5,4) = 1;
out(:,21,5,5) = 2;
end