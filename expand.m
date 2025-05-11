function [FF] = expand(F)
%EXPAND 此处显示有关此函数的摘要
%   此处显示详细说明
[n,~]=size(F);
if n==1
    F=F';
end

[n,~]=size(F);
c=length(unique(F));

FF=zeros(n,c);
for i=1:n
   FF(i,F(i))=1; 
end


end

