function [F,Flabel] = updataF(Flabel,A)
    [n,~] = size(Flabel);
    F = expand(Flabel);
    G = A + A';
    fTAf = diag(F'*A*F);
    fTf = diag(F'*F);
    fTg = F'*G;
    
    tol = 1e-6;
    maxcnt = 100;
    iterobjs = zeros(maxcnt,1);
    for iter = 1:maxcnt
        for i = 1:n
            r = Flabel(i);
            if fTf(r) == 1
                continue;
            end
            delta = (fTAf + fTg(:,i) + A(i,i)) ./(fTf + 1) - fTAf ./ fTf;
            delta(r) = fTAf(r) / fTf(r) - (fTAf(r) - fTg(r,i) + A(i,i)) /(fTf(r) - 1);
            [value, p]=max(delta);
            if p ~= r
                fTAf(r) = fTAf(r) - fTg(r,i) + A(i,i);
                fTAf(p) = fTAf(p) + fTg(p,i) + A(i,i);
                fTg(r,:) = fTg(r,:) - G(i,:);
                fTg(p,:) = fTg(p,:) + G(i,:);
                fTf(r) = fTf(r) - 1;
                fTf(p) = fTf(p) + 1;
                F(i,r) = 0;
                F(i,p) = 1;
                Flabel(i) = p;
            end
        end
        iterobjs(iter) = trace(inv(F'*F)*F'*A*F);      
        if (iter > 1) && abs((iterobjs(iter)-iterobjs(iter-1))/iterobjs(iter-1)) < tol
            break;
        end
    end
end

