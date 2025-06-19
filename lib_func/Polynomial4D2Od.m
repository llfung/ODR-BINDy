function out = Polynomial4D2Od(X)
N = numel(X)/4;
X = reshape(X,N,4);
   x1=X(:,1);
   x2=X(:,2);
   x3=X(:,3);
   x4=X(:,4);
out = zeros(N,15,4);
out(:,:,1) = [zeros(N,1) ones(N,1) zeros(N,1) zeros(N,1) zeros(N,1) 2.*x1 x2 x3 x4 zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) ];
out(:,:,2) = [zeros(N,1) zeros(N,1) ones(N,1) zeros(N,1) zeros(N,1) zeros(N,1) x1 zeros(N,1) zeros(N,1) 2.*x2 x3 x4 zeros(N,1) zeros(N,1) zeros(N,1) ];
out(:,:,3) = [zeros(N,1) zeros(N,1) zeros(N,1) ones(N,1) zeros(N,1) zeros(N,1) zeros(N,1) x1 zeros(N,1) zeros(N,1) x2 zeros(N,1) 2.*x3 x4 zeros(N,1) ];
out(:,:,4) = [zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) ones(N,1) zeros(N,1) zeros(N,1) zeros(N,1) x1 zeros(N,1) zeros(N,1) x2 zeros(N,1) x3 2.*x4 ];
end