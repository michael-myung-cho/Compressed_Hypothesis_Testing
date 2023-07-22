function [ BInd ] = HT_MP(Aseed, Y, n, ks, nc, mu, sigma, KList)
%% Compressed hypothesis testing Message Passing (MP) Ver3_1
% n: number of observation
% ni: seed matrix row dimension
% N: seed matrix col dimension
% ks: number of sparsity
% muO: mean of Odd (abnormal) RV
% muN: mean of Normal RV
% sigmaO: standard deviation of Odd (abnormal) RV
% sigmaN: standard deviation of Normal RV
% mu: mean of RV (mu = [muO,muN])
% sigma: standard deviation of RV (sigma = [sigmaO, sigmaN])

%% MP
[ni,N] = size(Aseed);
% yrMax = n*(ni >= n) + ni*(ni < n);
ycMax = ceil(n/ni);
muO = mu(1);
muN = mu(2);
sigmaO = sigma(1);
sigmaN = sigma(2);

pf = zeros(2^nc,ycMax,ni);
for j = 1:ni
    c1 = find(Aseed(j,:) ~= 0);
    ctA0 = [];
    ctA1 = [];
    for ii = 1:length(c1)
        ctA0 = [1*ones(2^(ii-1),1),ctA0; 0*ones(2^(ii-1),1),ctA0];
        ctA1 = [0*ones(2^(ii-1),1),ctA1; 1*ones(2^(ii-1),1),ctA1];
    end
    for l = 1:ycMax
        for ii = 1:2^length(c1)
            n1 = sum(ctA1(ii,:));
            pf(ii,l,j) = normpdf(Y(j,l), muO*n1 + muN*(length(c1)-n1), sqrt(sigmaO^2*n1 + sigmaN^2*(length(c1)-n1)));
        end
    end
end
% Initialization
[ qij0, qij1 ] = normpdfA2( Y, Aseed, sigma, mu, N, pf );
rji0 = zeros(ni,N);
rji1 = zeros(ni,N);

% Iteration
MaxItr = 20;
for itr = 1:MaxItr
    % fprintf('Iteration : %d\n', itr);
    % ----- Horizontal step -----
    for i = 1:N
        % Find non-zeros in the column
        r1 = find(Aseed(:, i) ~= 0);
        for k = 1:length(r1)
            drji0 = 1;
            drji1 = 1;
            for l=1:length(r1)
                if l ~= k
                    drji0 = drji0*qij0(r1(l),i);
                    drji1 = drji1*qij1(r1(l),i);
                end
            end
            Q0 = drji0/(drji1 + drji0);
            Q1 = drji1/(drji0 + drji1);
            
            rji0(r1(k),i) = Q0;
            rji1(r1(k),i) = Q1;
            
        end % for k
        % Update constants
        Qi0(i) = prod(qij0(r1,i));
        Qi1(i) = prod(qij1(r1,i));
        
        % Decode Qj
        if Qi1(i) > Qi0(i)
            vHat(i) = 1;
        else
            vHat(i) = 0;
        end
    end % for i
    
    
    % ------ Vertical step ------
    K0 = zeros(ni,N);
    K1 = zeros(ni,N);
    for j = 1:ni
        % Find non-zeros in the row
        c1 = find(Aseed(j,:) ~= 0);
        lc1 = length(c1);
        Ri1 = ctA1.*repmat(rji1(j, c1),2^lc1,1);
        Ri0 = ctA0.*repmat(rji0(j, c1),2^lc1,1);
        Ri = Ri1 + Ri0;
        nci = sum( Y(j,:) < 10^10 );
        for k = 1:lc1
            Rif = Ri;
            Rif(:,k) = 1;
            Rift = [Rif,pf(:,:,j)];
            pA = prod(Rift(:,1:(lc1+nci)),2); % row element multiplication
            
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
        end % for j
    end % for n
end
[Q1p,Q1ind] = sort(Qi1,'descend');
BInd = Q1ind(1:ks);

end




