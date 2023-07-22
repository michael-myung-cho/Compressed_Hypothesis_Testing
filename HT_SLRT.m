function [BInd] = HT_SLRT(XArray, n, N, ks, mu, sigma, KList)
%% Separate hypothesis testing Likelihood Ratio Test (SLRT)
% n: number of observation
% ni: seed matrix row dimension
% N: seed matrix col dimension
% ks: number of sparsity
% E1: mean of normal RV
% E2: mean of abnormal RV
% V1: variance of normal RV
% V2: variance of abnormal RV

nCase = nchoosek(N,ks);
E1 = mu(2); % normal
E2 = mu(1); % abnormal
V1 = sigma(2).^2; % normal
V2 = sigma(1).^2; % abnormal

%%% separate observations
%% generate separate observations
SepY=zeros(n,1);
SepA=zeros(n,N);
for i=1:n
    SepA(i,i-(ceil(i/N)-1)*N )=1;
end
%
%         for j=1:k
%             l=1;
%             while n*(l-1)+KList(1, k)<=m
%                 SepY(n*(l-1)+KList(1, k),1)=sqrt(V2)*randn+E2;
%                 l=l+1;
%             end
%         end
for i=1:n
    SepY(i,1)=SepA(i,:)*XArray(:,i);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%likelihood test
SepSet=nchoosek([1:N], ks);
Seplikelihood=zeros(nCase, 1);

for l=1:nCase
    
    
    SepYE=zeros(n,1);
    SepYV=zeros(n,1);
    %%%%%%%%%%%figure out the corresponding distribution
    %%%figure out the X distribution
    SepActiveSet=SepSet(l,:);
    SepEX=E1*ones(N,1);
    for j=1:ks
        SepEX(SepActiveSet(1,j),:)=E2;
    end
    
    SepVX=V1*ones(N,1);
    for j=1:ks
        SepVX(SepActiveSet(1,j),:)=V2;
    end
    
    
    %%%%%%%%%%%% Figure out the Y distirbution
    for p=1:n
        SepYE(p,1)=SepA(p,:)*SepEX;
        SepYV(p,1)=SepA(p,:).^2*SepVX;
    end
    
    %%%%calculate the likelihood
    
    for p=1:n
        Seplikelihood(l, 1)=Seplikelihood(l, 1)-0.5*log(2*pi*SepYV(p,1))-(SepY(p,1)-SepYE(p,1))^2/(2*SepYV(p,1));
    end
    
    
end

[Sepbestvalue, Sepbestindex]=max(Seplikelihood);
BInd=SepSet(Sepbestindex,:);
% if sort(SepSet(Sepbestindex,:))~=sort(KList)
%     SepEr(1,countm)=SepEr(1,countm)+1;
% end

end


