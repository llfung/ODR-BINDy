function out = Polynomial2D3Odd(X)
    N = numel(X)/2;
    X = reshape(X,N,2);
    out = zeros(N,10,2,2);

    % 2nd order terms
    out(:,4,1,1) = 2;
    out(:,5,2,1) = 1;
    out(:,5,1,2) = 1;
    out(:,6,2,2) = 2;
    % 3rd order terms
    out(:,7,1,1) = 6*X(:,1);
    out(:,8,1,2) = 2*X(:,1);
    out(:,8,2,1) = 2*X(:,1);
    out(:,9,1,2) = 2*X(:,2);
    out(:,9,2,1) = 2*X(:,2);
    out(:,10,2,2) = 6*X(:,2);

end