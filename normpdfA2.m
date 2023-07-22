function [ qij0, qij1 ] = normpdfA2( y, Aseed, sigma, mu, N, pf )
%NORMDP Summary of this function goes here
% Abnormal node mean and variable (m, sigma)
% normal node mean and variable (0, 1)
% n: number of nodes
sigmaO = sigma(1);
sigmaN = sigma(2);
muO = mu(1);
muN = mu(2);
[ni, nc] = size(y);
K0 = zeros(ni,N);
K1 = zeros(ni,N);
c1 = find(Aseed(1,:) ~= 0);
Ri = ones(2^length(c1),length(c1))/2;
ctA0 = [];
ctA1 = [];
for ii = 1:length(c1)
    ctA0 = [1*ones(2^(ii-1),1),ctA0; 0*ones(2^(ii-1),1),ctA0];
    ctA1 = [0*ones(2^(ii-1),1),ctA1; 1*ones(2^(ii-1),1),ctA1];
end
for j = 1:ni
    
    % Find non-zeros in the row
    c1 = find(Aseed(j,:) ~= 0);
    nci = sum( y(j,:) < 10^10 );
    for k = 1:length(c1)
        Rif = Ri;
        Rif(:,k) = 1;
        Rift = [Rif,pf(:,:,j)];
        pA = prod(Rift(:,1:length(c1)+nci),2);
        
        % Get row products of prodOfrij\ri(l)
        cs1 = find(ctA1(:,k) == 1);
        cs0 = find(ctA1(:,k) == 0);
        K0s = sum(pA(cs0));
        K1s = sum(pA(cs1));
        K0(j, c1(k)) = K0s;
        K1(j, c1(k)) = K1s;
        
        % Update qij0 and qij1
        qij0(j, c1(k)) = K0(j, c1(k))/(K1(j, c1(k)) + K0(j, c1(k)));
        qij1(j, c1(k)) = K1(j, c1(k))/(K0(j, c1(k)) + K1(j, c1(k)));
    end % for k
end % for j

end

