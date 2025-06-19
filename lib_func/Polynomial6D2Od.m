function out = Polynomial6D2Od(X)
N = numel(X)/6;
X = reshape(X,N,6);
   x1=X(:,1);
   x2=X(:,2);
   x3=X(:,3);
   x4=X(:,4);
   x5=X(:,5);
   x6=X(:,6);
out = zeros(N,28,6);
out(:,:,1) = [zeros(N,1) ones(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) 2.*x1 x2 x3 x4 x5 x6 zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) ];
out(:,:,2) = [zeros(N,1) zeros(N,1) ones(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) x1 zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) 2.*x2 x3 x4 x5 x6 zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) ];
out(:,:,3) = [zeros(N,1) zeros(N,1) zeros(N,1) ones(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) x1 zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) x2 zeros(N,1) zeros(N,1) zeros(N,1) 2.*x3 x4 x5 x6 zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) ];
out(:,:,4) = [zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) ones(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) x1 zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) x2 zeros(N,1) zeros(N,1) zeros(N,1) x3 zeros(N,1) zeros(N,1) 2.*x4 x5 x6 zeros(N,1) zeros(N,1) zeros(N,1) ];
out(:,:,5) = [zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) ones(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) x1 zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) x2 zeros(N,1) zeros(N,1) zeros(N,1) x3 zeros(N,1) zeros(N,1) x4 zeros(N,1) 2.*x5 x6 zeros(N,1) ];
out(:,:,6) = [zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) ones(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) x1 zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) x2 zeros(N,1) zeros(N,1) zeros(N,1) x3 zeros(N,1) zeros(N,1) x4 zeros(N,1) x5 2.*x6 ];
end