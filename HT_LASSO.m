function [ espInd ] = HT_LASSO(A, Y, ks, mu, sigma, KList)
%% Compressed hypothesis testing LASSO type 
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

[n,N]=size(A);
muO = mu(1);
muN = mu(2);
sigmaO = sigma(1);
sigmaN = sigma(2);

%% CVX implementation
lambda = 0.5;
c = zeros(N,1);

cvx_quiet true

cvx_begin
    variable c(N);
    minimize (sum_square(Y - A*c) + lambda*norm(c,1))
cvx_end
%% LASSO in Matlab implementation
% c = lasso(A,y,'Lambda',0.5);

[~,espIndA] = sort(abs(c),'descend');  % case when only the k biggest elements are considered
espInd = sort(espIndA(1:ks)');
end 