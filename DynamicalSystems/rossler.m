function dy = rossler(t,y,param,b,c)
    if nargin>3
        a=param;
    else
        a=param.a;
        b=param.b;
        c=param.c;
    end
    
    dy = [
    -y(3)-y(2);
    y(1)+a*y(2);
    b+y(3)*(y(1)-c);
    ];
end