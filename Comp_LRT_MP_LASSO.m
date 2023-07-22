%% Comparison: Likelihood Ratio Test (LRT) and Message Passing (MP) method, and LASSO type method
%
% (a) HT_CLRT: Hypothesis Testing via Compressed Likelihood Ratio Test 
% (b) HT_SLRT: Hypothesis Testing via Standard (Singe) Likelihood Ratio Test
% (c) HT_MP: Hypothesis Testing via Message Passing method
% (d) HT_LASSO: Hypothesis Testing via LASSO type method
%
% Sep. 23, 2016
% Myung (Michael) Cho (myung-cho@uiowa.edu)


clear
clc

nr = 8; % number of 1's in row of LDPC matrix
nc = 4; % number of 1's in col of LDPC matrix
ni = 200;
Aseed = genH_regularGallagher(ni, nr, nc);
% % genH_regularGallagher(n, row, col)
% % 1. (n/row)>col
% % 2. col>=3
% % 3. row>col
% % 4. row mod n must be zero

[ni,N] = size(Aseed)

%% save data to xls file
countXLS=1;    
clk=clock;
datestmp=strcat([date,'-',num2str(clk(4)),'h',num2str(clk(5))]);
filename = strcat('(',num2str(ni),'x',num2str(N),')','HT-Comp_LRT_MP_LASSO',datestmp,'_mean_0_8_trial1000.xlsx');
col_name = {'ni','N','k','n','pErrCLRT','pErrSLRT','pErrMP','pErrLASSO','nErrCLRT','nErrSLRT','nErrMP','nErrLASSO';};
xlswrite(filename,col_name);

%% two different variable definition
% mu: mean, sigma: standard deviation
muO = 8;
muN = 0;
mu = [muO; muN];
sigmaO = 10;
sigmaN = 1;
sigma = [sigmaO; sigmaN];

%% observation and recovery
nTrial = 1000;
ErrCLRTA = [];
ErrSLRTA = [];
ErrMPA = [];
iN = [];
nObse = 200;
kF = 1:1:2;
for ks = kF
    for n = ni:10:nObse
        fprintf('nObservation : %d\n', n);
        ErrCLRT = 0;
        ErrSLRT = 0;
        ErrMP = 0;
        ErrLASSO = 0;
        for iTrial = 1:nTrial
            fprintf('iTrial : %d\n', iTrial);
            %% observation y and deterministic sparse signal x
            Ind = randperm(N);
            KList = Ind(1:ks);
            x = [];
            A = [];
            Y = [];
            yrMax = n*(ni >= n) + ni*(ni < n);
            ycMax = ceil(n/ni);
            y = 10^10*ones(yrMax, ycMax); % initialization 
            for ii = 1:n
                xprime = muN + randn(N,1);
                xprime(Ind(1:ks)) = muO + sigmaO*randn(ks,1);
                yr = mod(ii-1,ni)+1;
                yc = ceil(ii/ni);
                x=[x,xprime];
                y(yr,yc) = Aseed(yr,:)*xprime;
                A = [A;Aseed(yr,:)];
                Y(ii,1) = y(yr,yc);
            end
            
            %% algorithm 
            [BIndCLRT] = HT_CLRT(A, Y, ks, mu, sigma, KList);
            if sum(sort(BIndCLRT)~=sort(KList))
                ErrCLRT=ErrCLRT+1;
            end
            [BIndSLRT] = HT_SLRT(x, n, N, ks, mu, sigma, KList );
            if sum(sort(BIndSLRT)~=sort(KList))
                ErrSLRT=ErrSLRT+1;
            end
            [BIndMP] = HT_MP(Aseed, y, n, ks, nc, mu, sigma, KList);
            if sum(sort(BIndMP)~=sort(KList))
                ErrMP=ErrMP+1;
            end
            [BIndLASSO] = HT_LASSO(A, Y, ks, mu, sigma, KList);
            if sum(sort(BIndLASSO)~=sort(KList))
                ErrLASSO=ErrLASSO+1;
            end
        end
        
        %% save data to xls file
        resBuf = [ni, N, ks, n, ErrCLRT/nTrial, ErrSLRT/nTrial, ErrMP/nTrial, ErrLASSO/nTrial, ErrCLRT, ErrSLRT, ErrMP, ErrLASSO];
        countXLS = countXLS+1;
        xlRange = sprintf('A%d',countXLS);
        xlsInsert = resBuf;
        xlswrite(filename, xlsInsert, 'Sheet1', xlRange);
    
    end
end

