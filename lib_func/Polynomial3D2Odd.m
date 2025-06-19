function out = Polynomial3D2Odd(X)
    N = numel(X)/3;
    % X = reshape(X,N,3);
    out = zeros(N,10,3,3);
    out(:,5,1,1) = 2;
    out(:,6,2,1) = 1;
    out(:,7,3,1) = 1;
    out(:,6,1,2) = 1;
    out(:,8,2,2) = 2;
    out(:,9,3,2) = 1;
    out(:,7,1,3) = 1;
    out(:,9,2,3) = 1;
    out(:,10,3,3) = 2;
end