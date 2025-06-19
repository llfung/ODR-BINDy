function out = Polynomial3D5Od(X)
N = numel(X)/3;
X = reshape(X,N,3);
   x1=X(:,1);
   x2=X(:,2);
   x3=X(:,3);
out = zeros(N,56,3);
out(:,:,1) = [zeros(N,1) ones(N,1) zeros(N,1) zeros(N,1) 2.*x1 x2 x3 zeros(N,1) zeros(N,1) zeros(N,1) 3.*x1.^2 2.*x1.*x2 2.*x1.*x3 x2.^2 x2.*x3 x3.^2 zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) 4.*x1.^3 3.*x1.^2.*x2 3.*x1.^2.*x3 2.*x1.*x2.^2 2.*x1.*x2.*x3 2.*x1.*x3.^2 x2.^3 x2.^2.*x3 x2.*x3.^2 x3.^3 zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) 5.*x1.^4 4.*x1.^3.*x2 4.*x1.^3.*x3 3.*x1.^2.*x2.^2 3.*x1.^2.*x2.*x3 3.*x1.^2.*x3.^2 2.*x1.*x2.^3 2.*x1.*x2.^2.*x3 2.*x1.*x2.*x3.^2 2.*x1.*x3.^3 x2.^4 x2.^3.*x3 x2.^2.*x3.^2 x2.*x3.^3 x3.^4 zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) ];
out(:,:,2) = [zeros(N,1) zeros(N,1) ones(N,1) zeros(N,1) zeros(N,1) x1 zeros(N,1) 2.*x2 x3 zeros(N,1) zeros(N,1) x1.^2 zeros(N,1) 2.*x1.*x2 x1.*x3 zeros(N,1) 3.*x2.^2 2.*x2.*x3 x3.^2 zeros(N,1) zeros(N,1) x1.^3 zeros(N,1) 2.*x1.^2.*x2 x1.^2.*x3 zeros(N,1) 3.*x1.*x2.^2 2.*x1.*x2.*x3 x1.*x3.^2 zeros(N,1) 4.*x2.^3 3.*x2.^2.*x3 2.*x2.*x3.^2 x3.^3 zeros(N,1) zeros(N,1) x1.^4 zeros(N,1) 2.*x1.^3.*x2 x1.^3.*x3 zeros(N,1) 3.*x1.^2.*x2.^2 2.*x1.^2.*x2.*x3 x1.^2.*x3.^2 zeros(N,1) 4.*x1.*x2.^3 3.*x1.*x2.^2.*x3 2.*x1.*x2.*x3.^2 x1.*x3.^3 zeros(N,1) 5.*x2.^4 4.*x2.^3.*x3 3.*x2.^2.*x3.^2 2.*x2.*x3.^3 x3.^4 zeros(N,1) ];
out(:,:,3) = [zeros(N,1) zeros(N,1) zeros(N,1) ones(N,1) zeros(N,1) zeros(N,1) x1 zeros(N,1) x2 2.*x3 zeros(N,1) zeros(N,1) x1.^2 zeros(N,1) x1.*x2 2.*x1.*x3 zeros(N,1) x2.^2 2.*x2.*x3 3.*x3.^2 zeros(N,1) zeros(N,1) x1.^3 zeros(N,1) x1.^2.*x2 2.*x1.^2.*x3 zeros(N,1) x1.*x2.^2 2.*x1.*x2.*x3 3.*x1.*x3.^2 zeros(N,1) x2.^3 2.*x2.^2.*x3 3.*x2.*x3.^2 4.*x3.^3 zeros(N,1) zeros(N,1) x1.^4 zeros(N,1) x1.^3.*x2 2.*x1.^3.*x3 zeros(N,1) x1.^2.*x2.^2 2.*x1.^2.*x2.*x3 3.*x1.^2.*x3.^2 zeros(N,1) x1.*x2.^3 2.*x1.*x2.^2.*x3 3.*x1.*x2.*x3.^2 4.*x1.*x3.^3 zeros(N,1) x2.^4 2.*x2.^3.*x3 3.*x2.^2.*x3.^2 4.*x2.*x3.^3 5.*x3.^4 ];
end