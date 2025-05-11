clear all;


name = "binalpha"

filename=strcat("bc_pool_",name);
load(filename);
k=length(unique(gt));
M=20;

members = members(:,1:M);
[n,~]=size(members);
epsilon=3;

cs=zeros(M,1);  
for m=1:M
    cs(m,1)=length(unique(members(:,m)));
end

% compute microH
microCls_num=sum(cs);
microH=zeros(n,microCls_num);
mH_index=1;
for m=1:M
    microC_num=cs(m);
    temp=expand(members(:,m));
    microH(:,mH_index:mH_index+microC_num-1)=temp;
    mH_index=mH_index+microC_num;
end

%compute P
p=zeros(microCls_num,microCls_num);
for i=1:microCls_num
    for j=1:microCls_num
        p(i,j)=(microH(:,i)'*microH(:,j))^2;
    end
end

ROUND = 15;
results=zeros(ROUND,3);
for round = 1:ROUND
% inite w
w=ones(microCls_num,1)./microCls_num;

%inite A
A=zeros(n,n);
for i=1:microCls_num
    temp=w(i)*microH(:,i)*microH(:,i)';
    A=A+temp;
end

%inite F
Flabel=randi([1,k],n,1);
F=expand(Flabel);

oldOBJ = getOBJ(w,A,F,microH,epsilon);
maxNITER = 6;
for iter=1:maxNITER
    %updata F
    % [F,Flabel]=get_newF2(A,k,Flabel);
    [F,Flabel] = updataF(Flabel,A);

    %compute e
    e=zeros(microCls_num,1);
    for i=1:microCls_num
        e(i)=microH(:,i)'*A*microH(:,i)*2;
    end

    %updata w
    Aeq=ones(1,microCls_num);
    beq=1;
    lb=zeros(microCls_num,1);
    [w,fval,exitflag]=quadprog(p,e,[],[],Aeq,beq,lb);

    %updata A
    S=zeros(n,n);
    A=zeros(n,n);
    for i=1:microCls_num
        temp=w(i)*microH(:,i)*microH(:,i)';
        S=S+temp;
    end
    B=F*inv(F'*F)*F';
    V=(2*S+epsilon*B)/2;
    for i=1:n
        A(i,:)=EProjSimplex_new(V(i,:));
    end
    newOBJ = getOBJ(w,A,F,microH,epsilon);
    
    if abs(newOBJ-oldOBJ) < 1e-6
        break;
    end
    oldOBJ = newOBJ;
    res=ClusteringMeasure(Flabel,gt);
    results(round,:)=res;
end
end
st=std(results);
me=mean(results);
