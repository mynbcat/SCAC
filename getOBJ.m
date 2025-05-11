function [OBJ] = getOBJ(w,A,F,microH,epsilon)
    [microCls_num,~] = size(w);
    [n,~] = size(A);
    S = zeros(n,n);
    for i=1:microCls_num
        temp=w(i)*microH(:,i)*microH(:,i)';
        S=S+temp;
    end
    output=norm(A-S,'fro')^2 - epsilon*trace(inv(F'*F)*F'*A*F);
    OBJ = output;
end

