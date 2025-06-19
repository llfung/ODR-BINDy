function out = Polynomial2D3O(X)
    N = numel(X)/2;
    X = reshape(X,N,2);
    out = [ones(N,1) X(:,1) X(:,2) ...
        X(:,1).^2 X(:,1).*X(:,2) X(:,2).^2 ...
        X(:,1).^3 X(:,1).^2.*X(:,2) X(:,2).^2.*X(:,1) X(:,2).^3];
end