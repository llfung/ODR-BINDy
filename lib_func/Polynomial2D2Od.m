function out = Polynomial2D2Od(X)
N = numel(X)/2;
X = reshape(X,N,2);
   x1=X(:,1);
   x2=X(:,2);
out = zeros(N,6,2);
out(:,:,1) = [zeros(N,1) ones(N,1) zeros(N,1) 2.*x1 x2 zeros(N,1) ];
out(:,:,2) = [zeros(N,1) zeros(N,1) ones(N,1) zeros(N,1) x1 2.*x2 ];
end