function [ BInd ] = HT_CLRT(A, Y, ks, mu, sigma, KList)
%% Compressed hypothesis testing Likelihood Ratio Test (CLRT)
% n: number of observation
% ni: seed matrix row dimension
% N: seed matrix col dimension
% ks: number of sparsity
% E1: mean of normal RV
% E2: mean of abnormal RV
% V1: variance of normal RV
% V2: variance of abnormal RV

[n,N]=size(A);

E1=mu(2);  % normal 
E2=mu(1);  % abnormal
V1=sigma(2).^2;  % normal
V2=sigma(1).^2;  % abnormal

%% likelihood test
Set=nchoosek([1:N], ks);
nCase=nchoosek(N,ks);
likelihood=zeros(nCase, 1);

for l=1:nCase

    YE=zeros(n,1);
    YV=zeros(n,1);
    
    %%%%%%%%%%%figure out the corresponding distribution
    %%%figure out the X distribution
    ActiveSet=Set(l,:);
    EX=E1*ones(N,1);
    for j=1:ks
        EX(ActiveSet(1,j),:)=E2;
    end
    
    VX=V1*ones(N,1);
    for j=1:ks
        VX(ActiveSet(1,j),:)=V2;
    end
    
    
    %%%%%%%%%%%% Figure out the Y distirbution
    for p=1:n
        YE(p,1)=A(p,:)*EX;
        YV(p,1)=A(p,:).^2*VX;
    end
    
    %%%%calculate the likelihood
    
    for p=1:n
        likelihood(l, 1)=likelihood(l, 1)-0.5*log(2*pi*YV(p,1))-(Y(p,1)-YE(p,1))^2/(2*YV(p,1));
    end
    
    
end

[bestvalue, bestindex]=max(likelihood);
BInd=Set(bestindex,:);

end


