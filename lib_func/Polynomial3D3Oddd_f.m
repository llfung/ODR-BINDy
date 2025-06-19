
function out = Polynomial3D3Oddd_f(X,p_selected,mask_locd)
    % output = Theta'''(x)*p given only mask_locd term active
    N = numel(X)/3;
    out = zeros(N,3,3,3);
    p = zeros(20,1);
    p(mask_locd) = p_selected;

    out(:,1,1,1) = 6*p(11);
    out(:,1,1,2) = 2*p(12);
    out(:,1,2,1) = 2*p(12);
    out(:,2,1,1) = 2*p(12);
    out(:,1,1,3) = 2*p(13);
    out(:,1,3,1) = 2*p(13);
    out(:,3,1,1) = 2*p(13);
    out(:,1,2,2) = 2*p(14);
    out(:,2,1,2) = 2*p(14);
    out(:,2,2,1) = 2*p(14);
    out(:,1,2,3) = 2*p(15);
    out(:,3,1,2) = 2*p(15);
    out(:,2,3,1) = 2*p(15);
    out(:,3,2,1) = 2*p(15);
    out(:,2,1,3) = 2*p(15);
    out(:,1,3,2) = 2*p(15);
    out(:,1,3,3) = 2*p(16);
    out(:,3,1,3) = 2*p(16);
    out(:,3,3,1) = 2*p(16);
    out(:,2,2,2) = 6*p(17);
    out(:,2,2,3) = 2*p(18);
    out(:,2,3,2) = 2*p(18);
    out(:,3,2,2) = 2*p(18);
    out(:,2,3,3) = 2*p(19);
    out(:,3,2,3) = 2*p(19);
    out(:,3,3,2) = 2*p(19);   
    out(:,3,3,3) = 6*p(20);
end