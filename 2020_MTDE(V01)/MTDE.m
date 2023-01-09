%___________________________________________________________________________________%
%  MTDE: An effective multi-trial vector-based differential evolution algorithm     %
%  and its applications for engineering design problems                             %
%  source codes version 1.0                                                         %
%                                                                                   %
%  Developed in MATLAB R2018a                                                       %
%                                                                                   %
%  Author and programmer: M. H.Nadimi-Shahraki, S. Taghian, S. Mirjalili, H. Faris  %
%                                                                                   %
%  e-Mail:nadimi@ieee.org, shokooh.taghian94@gmail.com,                             %                                                                                          
%  ali.mirjalili@gmail.com, hossam.faris@ju.edu.jo                                  %
%                                                                                   %
%  Homepage:http://www.alimirjalili.com                                             %
%                                                                                   %
%  Main paper: M. H. Nadimi-Shahraki, S. Taghian, S. Mirjalili, H. Faris            %
%  MTDE: An effective multi-trial vector-based differential evolution               %
%  algorithm and its applications for engineering design problems,                  %
%  Applied Soft Computing Journal 97 (2020) 106761                                  %
%                                                                                   %
%  DOI: 10.1016/j.asoc.2020.106761                                                  %
%___________________________________________________________________________________%

function [gbestval,gbest,BestChart] = MTDE(D,N,MaxFES,lb,ub,fobj)
rand('seed', sum(100 * clock));

%% the values and indices of the best solutions
FES = 0;
Gen = 0;
WinIter = 20;
nFES = [0,0,0];
MaxGen = MaxFES/N;

winTVPIdx = 1;
nImpFit = [1,1,1];
ImpRate = [0,0,0];
leastPortion = 0.2;
sf_memory_size = D;
F = 0.5;

H = 5;                               % The history size of the gbest archive
lu = [lb * ones(1, D); ub * ones(1, D)];

BestChart=[];

%% Initialization
X = repmat(lu(1, :), N, 1) + rand(N, D) .* (repmat(lu(2, :) - lu(1, :), N, 1));
for i=1:N
    Fit(i) = fobj(X(i,:));
end

[gbestval,tmp1] = min(Fit);
gbest = X(tmp1, :);

gb_history_X = X(tmp1,:);            % gbest archive
gb_history_Fit = gbestval;
Mgb_h = repmat(gbest,N,1);           % Mgb_h::X_hat_gbest
count = 1;

posbest = X;
pbestval = Fit;

  %% ======================================================================%%%%
memory_sf = 0.5 .* ones(sf_memory_size, 1);
memory_pos = 1;
%
LifeTime_archive.N = N;              % The maximum size of the Life-time archive
LifeTime_archive.X = zeros(0, D);    % The stored individuals in the Life-time archive
LifeTime_archive.Fit = zeros(0, 1);  % The fitness value of the archived individuals
lifeTimeCounter = zeros(N, 1);
%
repeat = floor(N / D);
copy = mod(N, D);
%
initial = 0.001;
final = 2;
Mu = log(D);
%%
while  FES < MaxFES
    Gen = Gen +1;
    ca1 = 2 - Gen * ((2) /MaxGen);                                          % Eq. (9)
    ca2 = (initial - (initial - final) * (((MaxGen - Gen)/MaxGen))^Mu);     % Eq. (11)
    
    %% ===========================Determine scale factor (F) =====================================%%%%
    mem_rand_index = ceil(sf_memory_size * rand(N, 1));
    mu_sf = memory_sf(mem_rand_index);
    
    % Generating scaling factor
    F = mu_sf + 0.2 * tan(pi * (rand(N, 1) - 0.5));
    pos = find(F <= 0);
    
    while ~ isempty(pos)
        F(pos) = mu_sf(pos) + 0.2 * tan(pi * (rand(length(pos), 1) - 0.5));
        pos = find(F <= 0);
    end
    F = min(F, 1);
    
    %% ===========================Transformation Matrix (M) =====================================%%%%
    M_tril = tril(ones(D,D));
    if (copy == 0)&&(N ~= D)
        M_tril2= repmat(M_tril, repeat, 1);
    else
        M_tril2 = repmat(M_tril, repeat, 1);
        added_row = M_tril(1:copy,:);
        M_tril2 = [M_tril2;added_row];
    end
    
    Tmp = zeros(N,D);
    
    [nRows,nCols] = size(M_tril2);
    [~,idx] = sort(rand(nRows,nCols),2);
    idx = (idx - 1) * nRows + ndgrid(1:nRows,1:nCols);
    Tmp(:) = M_tril2(idx);
    
    M = Tmp(randperm(size(Tmp, 1)), :);
    M_bar = ~M;
    
    %%  %% ===========================Create Mgb_h =====================================%%%%
    nCopy_gbest = mod(N,count);
    rep = floor(N/count);
    
    if nCopy_gbest == 0
        Mgb_h = repmat(gb_history_X, rep, 1);
    else
        Mgb_h = repmat(gb_history_X, rep, 1);
        add = gb_history_X(1:nCopy_gbest,:);
        Mgb_h = [Mgb_h;add];
    end
    
    %%  %% ===========================The winner-based distributing substep=====================================%%%%
    if mod(Gen,WinIter) == 0                        % Definition 1/Eq. (2)
        ImpRate(1) = nImpFit(1)/nFES(1);
        ImpRate(2) = nImpFit(2)/nFES(2);
        ImpRate(3) = nImpFit(3)/nFES(3);
        
        [~,winTVPIdx] = max(ImpRate);
        
        nImpFit = [1,1,1];
        ImpRate =  [0,0,0];
        nFES = [0,0,0];
    end
    
    % Generate Subpopulations
    permutation = randperm(N);                      % Definition 2/Eqs. (3-4)
    if winTVPIdx == 1
        arrayGTVP= permutation(1:leastPortion*N);
        arrayLTVP = permutation(leastPortion*N+1:2*leastPortion*N);
        arrayRTVP = permutation(2*leastPortion*N+1:end);
    elseif winTVPIdx == 2
        arrayGTVP = permutation(1:leastPortion*N);
        arrayRTVP = permutation(leastPortion*N+1:2*leastPortion*N);
        arrayLTVP  = permutation(2*leastPortion*N+1:end);
    elseif winTVPIdx == 3
        arrayRTVP = permutation(1:2*leastPortion*N);
        arrayLTVP = permutation(2*leastPortion*N+1: 4*leastPortion*N);
        arrayGTVP  = permutation(4*leastPortion*N+1:end);
    end
    nFES = nFES + [length(arrayRTVP),length(arrayLTVP),length(arrayGTVP)];
    
    %%  %% ===========================R-TVP=====================================%%%%
    if ~isempty(arrayRTVP)
        X_RTVP = X(arrayRTVP,:);                    % The dedicated portion of the population to R-TVP
        Fit_RTVP = Fit(arrayRTVP);
        N_RTVP = length(arrayRTVP);
        F1 = F(arrayRTVP);
        
        r0 = 1 : N_RTVP;
        X_Upop = [X;LifeTime_archive.X];
        [r1, r2] = gnR1R2(N_RTVP, size(X_Upop, 1), r0);
        % Find the p-best solutions
        [~, indBest] = sort(Fit_RTVP, 'ascend');
        randindex = randi(N_RTVP,1,N_RTVP);
        pbest = X_RTVP(indBest(randindex), :);      % randomly choose one of the top 100p% solution
        
        % Find the p-worst solutions
        [~, indWorst] = sort(Fit_RTVP, 'descend');
        randindex = randi(N_RTVP,1,N_RTVP);
        pworst = X_RTVP(indWorst(randindex), :);    % randomly choose one of the top 100p% solution
        % == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==
        vi = X_RTVP + F1(:,ones(1,D)) .* (pbest - X_RTVP) + F1(:,ones(1,D)) .* (pworst - X_RTVP) + ...
            ca1 .* (X_Upop(r2,:) - X_RTVP);
        
        ui = M(arrayRTVP,:) .* X_RTVP + M_bar(arrayRTVP,:) .* vi;
        ui = boundConstraint(ui, lu);
        % == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==       
        for i=1:N_RTVP
            tmp_Fit_RTVP(i) = fobj(ui(i,:));
        end
        tmp = (Fit_RTVP <= tmp_Fit_RTVP);
        nImpFit(1) = nImpFit(1) + sum(tmp == 0);
        
        X_trial(arrayRTVP,:) = ui;
        Fit_X_trial(arrayRTVP) = tmp_Fit_RTVP;
    end
    %% ===========================L-TVP=====================================%%%%
    if ~isempty(arrayLTVP)
        X_LTVP = X(arrayLTVP,:);                    % The dedicated portion of the population to L-TVP
        Fit_LTVP = Fit(arrayLTVP);
        N_LTVP = length(arrayLTVP);
        F2 = F(arrayLTVP);
        
        r0 = 1 : N_LTVP;
        X_Upop = [X;LifeTime_archive.X];
        [r1, r2] = gnR1R2(N_LTVP, size(X_Upop, 1), r0);
        
        rot = (0:1:N_LTVP-1);                       % Shuffled population for pm1 and pm2
        ind = randperm(2);
        a1  = randperm(N_LTVP);                     
        rt = rem(rot+ind(1),N_LTVP);                
        a2  = a1(rt+1);                             
        rt = rem(rot+ind(2),N_LTVP);
        a3  = a2(rt+1);
        pm1 = X_LTVP(a2,:);                         
        pm2 = X_LTVP(a3,:);                         
        % == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==
        vi = X_LTVP + ca2.*(X_Upop(r2,:) - X_LTVP)+ F2(:,ones(1,D)).* (pm1 - pm2);
        
        ui = vi;
        ui = boundConstraint(ui, lu);
        % == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==
        for i=1:N_LTVP
            tmp_Fit_LTVP(i) = fobj(ui(i,:));
        end
        tmp = (Fit_LTVP <= tmp_Fit_LTVP);
        nImpFit(2) = nImpFit(2) + sum(tmp == 0);
        
        X_trial(arrayLTVP,:) = ui;
        Fit_X_trial(arrayLTVP) = tmp_Fit_LTVP;
    end
    
    %% ===========================G-TVP=====================================%%%%
    if ~isempty(arrayGTVP)
        X_GTVP = X(arrayGTVP,:);                    % The dedicated portion of the population to G-TVP
        Fit_GTVP = Fit(arrayGTVP);
        N_GTVP = length(arrayGTVP);
        
        r0 = 1 : N_GTVP;
        [r1, r2] = gnR1R2(N_GTVP, size(X_GTVP, 1), r0);
        X_differ = X_GTVP(r1, :) - X_GTVP(r2, :);
        
        vi = Mgb_h(arrayGTVP,:) + ca2.*  X_differ;
        
        ui = M(arrayGTVP,:) .* X_GTVP + M_bar(arrayGTVP,:) .* vi;
        ui = boundConstraint(ui, lu);
        
        % == == == == == == == == == == == == == == == == == == == == == == == == == == == ==        
        for i=1:N_GTVP
            tmp_Fit_GTVP(i) = fobj(ui(i,:));
        end
        
        tmp = (Fit_GTVP <= tmp_Fit_GTVP);
        nImpFit(3) = nImpFit(3) + sum(tmp == 0);
        
        X_trial(arrayGTVP,:) = ui;
        Fit_X_trial(arrayGTVP) = tmp_Fit_GTVP;
    end
    
    %% ======================= Population updating =========================================%%%%
    bin = (Fit_X_trial < pbestval)';
    diff = abs(pbestval - Fit_X_trial);
    
    goodF = F(bin == 1);
    dif_val = diff(bin == 1);
    
    [LifeTime_archive, lifeTimeCounter] = updateArchive(LifeTime_archive, posbest(bin == 1, :), pbestval(bin == 1)', lifeTimeCounter);
    
    posbest(bin==1,:) = X_trial(bin==1,:);
    pbestval(bin==1) = Fit_X_trial(bin==1);
    
    Sf = F(bin==1);
    diffB = diff(bin==1);
    numSucc = numel(Sf);
    %% ================================================================%%%%
    if numSucc > 0
        sum_dif = sum(dif_val);
        dif_val = dif_val / sum_dif;
        % Updating the memory of scaling factor
        memory_sf(memory_pos) = (dif_val * (goodF .^ 2)) / (dif_val * goodF);
        memory_pos = memory_pos + 1;
        if memory_pos > sf_memory_size;  memory_pos = 1; end
    end
    
    %% ================================================================%%%%
    [gbestval, gbestid] = min(pbestval);
    gbest = posbest(gbestid,:);
    
    mean_gbest_history = mean(gb_history_Fit);
    member = ~any((ismembertol(gb_history_Fit, gbestval))== 1);
    if (gbestval < mean_gbest_history)&&(member)
        if count >= H                               % The gbest archive is full
            [~, Idx] = max(gb_history_Fit);
            gb_history_Fit(Idx) = gbestval;
            gb_history_X(Idx,:) = gbest;
        else
            count = count+1;
            gb_history_X = [gb_history_X;gbest];
            gb_history_Fit = [gb_history_Fit;gbestval];
        end
    end
    %% ================================================================%%%%
    X = posbest;
    Fit = pbestval;
    BestChart = [BestChart gbestval];
    %     disp(['iter = ', num2str(Gen),' ,BEST SOLUTION MTDE= ' num2str(gbestval)]);
    
    clear tmp_Fit_RTVP tmp_Fit_LTVP tmp_Fit_GTVP;
    FES = FES + N;
    
end
end







