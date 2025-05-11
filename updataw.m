function [outputArg1,outputArg2] = updataw(w,P,e)
    [n,~] = size(w);
    mu = 0.5;
    rou = 1.5;
    multi = eye(n);
    old_objvalue = w'*P*w - w'*e;
    
    for i = 1:30
        beta = w + 1/mu*(multi - P'*w);
        r = 1/mu*(mu*beta - multi - P * beta + e);
        w = EProjSimplex_new(r);
        multi = multi + mu*(w - beta);
        mu = rou * mu;
        new_objvalue = w'*P*w - w'*e;
        if(abs((new_objvalue - old_objvalue)/old_objvalue) < e-4)
            break;
        end
        old_objvalue = new_objvalue;
    end
end

