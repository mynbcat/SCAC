function [newF,newLable] = CDEC(A,F,label,c)
%CDEC1 Summary of this function goes here
%   Detailed explanation goes here

[n,~] = size(A);
aa = sum(F,1);
last = 0;
iter_num = 0;


% caculate initial obj value
% for ii=1:c
%     if aa(ii) == 0
%         continue
%     end
%     
%     fl = F(:,ii);
%     temp = fl' * A * fl;
%     temp1 = aa(ii);
%     sumd(ii,1) = temp / temp1;
%     
% end
% OBJ(1) = sum(sumd);

OBJ(1) = getOBJ(A,F);

AF = A * F;
FAF = F'*A*F;
FTA = F' * A;

while any(label ~= last) && iter_num < 100
 last = label;       
 for i = 1:n   
    m = label(i) ;
    if aa(m)==1
        continue;  
    end 
    for k = 1:c        
        if k == m   
           V1(k) = FAF(k,k)-FTA(k,i)-AF(i,k)+A(i,i);
           delta(k) = FAF(k,k) / aa(k) - V1(k) / (aa(k) -1); 
        else  
           V2(k) = FAF(k,k)  + FTA(k,i)+AF(i,k)+A(i,i);
           delta(k) = V2(k) / (aa(k) +1) -  FAF(k,k)  / aa(k); 
        end         
    end  
    [~,q] = max(delta);     
    if m~=q
%           AF(:,q) = AF(:,q)+A(:,i);
%           AF(:,m) = AF(:,m)-A(:,i); 
          F(i,q) = 1;
          F(i,m) = 0;
          aa(q)= aa(q) + 1; %  FF(p,p)=F(:,p)'*F(:,p);
          aa(m)= aa(m) - 1; %  FF(m,m)=F(:,m)'*F(:,m)
%           FAF(m,m)=V1(m);
%           FAF(q,q)=V2(q);
%           FTA(m,i) = FTA(m,i)-A(m,i);
%           FTA(q,i) = FTA(q,i)+A(m,i);
          AF = A * F;
          FAF = F'*A*F;
          FTA = F' * A;
          label(i)=q;
    end
 end 
%     F = zeros(n,c);
%     for iiT = 1:n
%         F(iiT,label(iiT)) = 1;
%     end
    iter_num = iter_num + 1;
%     aa = sum(F,1);
%     for iiC=1:c
%         if aa(iiC) == 0
%             continue
%         end
% 
%         fl = F(:,iiC);
%         temp = fl' * A * fl;
%         temp1 = aa(iiC);
%         sumd(iiC,1) = temp / temp1;
%     
%     end
%     OBJ(iter_num+1) = sum(sumd);
    OBJ(iter_num+1) = getOBJ(A,F);
end
newLable = label;
newF = F;
end

