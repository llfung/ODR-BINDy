function out = Polynomial3D3O(X)
    N = numel(X)/3;
    X = reshape(X,N,3);
    out = [ones(N,1) X(:,1) X(:,2) X(:,3) ...
        X(:,1).^2 X(:,1).*X(:,2) X(:,1).*X(:,3) X(:,2).^2 X(:,2).*X(:,3) X(:,3).^2 ...
        X(:,1).^3 X(:,1).^2.*X(:,2) X(:,1).^2.*X(:,3) X(:,1).*X(:,2).^2 X(:,1).*X(:,2).*X(:,3) X(:,1).*X(:,3).^2 ...
        X(:,2).^3 X(:,2).^2.*X(:,3) X(:,2).*X(:,3).^2 X(:,3).^3];
end