
function out = Polynomial2D3Oddd_f(X,p_selected,mask_locd)
    % output = Theta'''(x)*p given only mask_locd term active
    N = numel(X)/2;
    out = zeros(N,2,2,2);
    p = zeros(10,1);
    p(mask_locd) = p_selected;

    out(:,1,1,1) = 6*p(7);
    out(:,2,1,1) = 2*p(8);
    out(:,1,2,1) = 2*p(8);
    out(:,1,1,2) = 2*p(8);
    out(:,2,2,1) = 2*p(9);
    out(:,1,2,2) = 2*p(9);
    out(:,2,1,2) = 2*p(9);
    out(:,2,2,2) = 6*p(10);
end