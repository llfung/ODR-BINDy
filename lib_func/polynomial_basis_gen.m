% Copyright 2025, All Rights Reserved
% Code by Lloyd Fung
%    For package ODR-BINDy
%

D = 10;
degree = 3;
syms x [1 D];

yout = [1.0];

for deg=1:degree
    yout = [yout polynomial_basis_at_deg(x,deg,D,1,[])];
end

f_name = strcat("Polynomial",num2str(D),"D",num2str(degree),"O");
df_name = strcat(f_name,"d");
ddf_name = strcat(f_name,"dd");
f2file(yout,f_name,D);
df2file(yout,df_name,D,x);
ddf2file(yout,ddf_name,D,x);
% matlabFunction(yout,"File",f_name,"Optimize",false);
% 
% fileID = fopen(strcat(f_name,".m"),'r+');
% topline = fgetl(fileID);
% while ~contains(topline,strcat("function yout = ",f_name))
%     topline = fgetl(fileID);
% end
% topline = 
%% Internal function
%Generate all polynomial terms at specific degree "deg"
% Recursive function: call with ii=1 and yout=[] in the beginning
function yout = polynomial_basis_at_deg(x,deg,D,ii,yout)
    if deg==1
        for i = ii:D
            yout = [yout x(i)];
        end
    else
        yout0=yout;
        for i = ii:D
            yout = [yout x(i).*polynomial_basis_at_deg(x,deg-1,D,i,yout0)];
        end
    end
end
% Generate Script with custom function (not using matlabFunction because of
% lack of support for vector input
function f2file(yout,f_name,D)
    fileID = fopen(strcat(f_name,".m"),'w');
    fprintf(fileID,strcat("function out = ",f_name,"(X)\n"));
    fprintf(fileID,strcat("N = numel(X)/",num2str(D),";\n"));
    fprintf(fileID,strcat("X = reshape(X,N,",num2str(D),");\n"));
    for i = 1:D
        fprintf(fileID,strcat("   x",num2str(i),"=X(:,",num2str(i),");\n"));
    end
    fprintf(fileID,"out = [");
    for i = 1:length(yout)
        if isSymType(yout(i),'number')
            if double(yout(i))==0
                fprintf(fileID,"zeros(N,1) ");
            elseif double(yout(i))==1
                fprintf(fileID,"ones(N,1) ");
            else 
                fprintf(fileID,strcat(num2str(double(yout(i))),"*ones(N,1) "));
            end
        else
            tmp = string(yout(i));
            tmp=replace(tmp,"^",".^");
            tmp=replace(tmp,"*",".*");
            tmp=replace(tmp,"/","./");
            fprintf(fileID,strcat(tmp," "));
        end
    end
    fprintf(fileID,"];\n");
    fprintf(fileID,"end");
    fclose(fileID);
end
function df2file(yout,f_name,D,x)
    fileID = fopen(strcat(f_name,".m"),'w');
    fprintf(fileID,strcat("function out = ",f_name,"(X)\n"));
    fprintf(fileID,strcat("N = numel(X)/",num2str(D),";\n"));
    fprintf(fileID,strcat("X = reshape(X,N,",num2str(D),");\n"));
    for i = 1:D
        fprintf(fileID,strcat("   x",num2str(i),"=X(:,",num2str(i),");\n"));
    end
    fprintf(fileID,strcat("out = zeros(N,",num2str(length(yout)),",",num2str(D),");\n"));

    for d = 1:D
        dyout = diff(yout,x(d));
        fprintf(fileID,strcat("out(:,:,",num2str(d),") = ["));
        for i = 1:length(yout)
            if isSymType(dyout(i),'number')
                if double(dyout(i))==0
                    fprintf(fileID,"zeros(N,1) ");
                elseif double(dyout(i))==1
                    fprintf(fileID,"ones(N,1) ");
                else 
                    fprintf(fileID,strcat(num2str(double(dyout(i))),"*ones(N,1) "));
                end
            else
                tmp = string(dyout(i));
                tmp=replace(tmp,"^",".^");
                tmp=replace(tmp,"*",".*");
                tmp=replace(tmp,"/","./");
                fprintf(fileID,strcat(tmp," "));
            end
        end
        fprintf(fileID,"];\n");
    end
    
    fprintf(fileID,"end");
    fclose(fileID);
end
function ddf2file(yout,f_name,D,x)
    fileID = fopen(strcat(f_name,".m"),'w');
    fprintf(fileID,strcat("function out = ",f_name,"(X)\n"));
    fprintf(fileID,strcat("N = numel(X)/",num2str(D),";\n"));
    fprintf(fileID,strcat("X = reshape(X,N,",num2str(D),");\n"));
    for i = 1:D
        fprintf(fileID,strcat("   x",num2str(i),"=X(:,",num2str(i),");\n"));
    end
    fprintf(fileID,strcat("out = zeros(N,",num2str(length(yout)),",",num2str(D),",",num2str(D),");\n"));

    for d1 = 1:D
        for d2 = 1:D
            dyout = diff(yout,x(d1),x(d2));
            for i = 1:length(yout)
                if isSymType(dyout(i),'number')
                    if double(dyout(i))==0
                        
                    else
                        fprintf(fileID,strcat("out(:,",num2str(i),",",num2str(d1),",",num2str(d2),") = ",...
                            num2str(double(dyout(i))),";\n"));
                    end
                else
                    tmp = string(dyout(i));
                    tmp=replace(tmp,"^",".^");
                    tmp=replace(tmp,"*",".*");
                    tmp=replace(tmp,"/","./");
                    fprintf(fileID,strcat("out(:,",num2str(i),",",num2str(d1),",",num2str(d2),") = ",...
                            tmp,";\n"));
                end
            end
        end
    end
    
    fprintf(fileID,"end");
    fclose(fileID);
end