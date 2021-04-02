%% The following code executes the LSADE algorithm
%% The algorithm is written in a maximization form (there are quite a few wild -1s appearing apparently out of nowhere)
%  The RBF surrogate model functions (rbf_build.m and rbf_predict.m) are from:
%       SURROGATES Toolbox: FAC Viana, \SURROGATES Toolbox User's Guide,"
%       Version 3.0, 2011, available at http://sites.google.com/site/felipeacviana/surrogatestoolbox.
%  The rest was created by: Jakub Kudela (Jakub.Kudela@vutbr.cz)

clc;clear;addpath('functions');
warning('off','all');

for ffff = [1] %choose the dynamic rule
    if ffff == 1
    % % 1-4 8-1
    f_li = @(iter) ceil((8*iter)/1000); f_lo = @(iter) ceil((8000-15*iter)/1000); sfx = ['_decr_Li14_Lo81'];
    elseif ffff == 2
    % 1-6 8-1
    f_li = @(iter) ceil((10*iter)/1000); f_lo = @(iter) ceil((8000-15*iter)/1000); sfx = ['_decr_Li16_Lo81'];
    elseif ffff == 3
    % % 1-8 8-1
    f_li = @(iter) ceil((14*iter)/1000); f_lo = @(iter) ceil((8000-15*iter)/1000); sfx = ['_decr_Li18_Lo81'];
    elseif ffff == 4    
    % % 1-4 6-1
    f_li = @(iter) ceil((8*iter)/1000); f_lo = @(iter) ceil((6000-12*iter)/1000); sfx = ['_decr_Li14_Lo61'];
    elseif ffff == 5
    % % 1-6 6-1
    f_li = @(iter) ceil((10*iter)/1000); f_lo = @(iter) ceil((6000-10*iter)/1000); sfx = ['_decr_Li16_Lo61'];
    elseif ffff == 6
    % % 1-8 6-1
    f_li = @(iter) ceil((14*iter)/1000); f_lo = @(iter) ceil((6000-10*iter)/1000); sfx = ['_decr_Li18_Lo61'];
    elseif ffff == 7
    % % 1-4 4-1
    f_li = @(iter) ceil((8*iter)/1000); f_lo = @(iter) ceil((4000-8*iter)/1000); sfx = ['_decr_Li14_Lo41'];
    elseif ffff == 8
    % % 1-6 4-1
    f_li = @(iter) ceil((12*iter)/1000); f_lo = @(iter) ceil((4000-8*iter)/1000); sfx = ['_decr_Li16_Lo41'];
    elseif ffff == 9
    % 1-8 4-1
    f_li = @(iter) ceil((15*iter)/1000); f_lo = @(iter) ceil((4000-8*iter)/1000); sfx = ['_decr_Li18_Lo41'];
    end

for D = [30]    % choice of dimensions
    for funcnr = [1:7]    % choice of function
        global initial_flag % needed to access the CEC 2005 functions
        initial_flag = 0;   % needed to access the CEC 2005 functions
        if funcnr == 1      % algorithm in maximization form, hence all functions need a '-'
           CostFunction = @(x) -sum((1:D).*x.^2); VarMin=-5.12; VarMax=5.12;
        elseif funcnr == 2
            CostFunction = @(x) -rosenbrock(x); VarMin=-2.048; VarMax=2.048;
        elseif funcnr == 3
            CostFunction = @(x) -ackley(x); VarMin=-32.768; VarMax=32.768;
        elseif funcnr == 4
            CostFunction = @(x) -griewank(x); VarMin=-600; VarMax=600;
        elseif funcnr == 5 
            CostFunction = @(x) -benchmark_func(x,10); VarMin=-5;VarMax= 5;
        elseif funcnr == 6 
            CostFunction = @(x) -benchmark_func(x,16); VarMin=-5;VarMax= 5;
        elseif funcnr == 7
            CostFunction = @(x) -benchmark_func(x,19); VarMin=-5;VarMax= 5;  
        end

        %% DE Parameters
        if D<=50
            initPop = 100;  % Initial Population Size
        else
            initPop = 200;
        end
        nPop = D;           % Population Size
        F=0.5;              % Scaling Factor
        CR=0.5;             % Crossover Probability
        %% Initialization
        c = 1; ker = 'MQ'; pol_part = 0;
        MaxIt = 1000;       % Maximum number of iterations (not function evaluations)
        MaxFe = 1000;       % Maximum number of function evaluations
        add_LIPO = 1; add_RBF = 1; add_FMINCON =1; % choice of adding individual surrogate based points
        best_results = zeros(1,20); best_vals = cell(1,20); nr_evals = cell(1,20); com_time = zeros(1,20);
        for trial = 1:20    % 20 runs
            tStart = tic;
            rng(trial,'twister');     % same initial conditions for comparison
            LIPO_flag = 0; RBF_flag = 0;
            X_all = (VarMax-VarMin)*lhsdesign(initPop,D)+VarMin; pop_cost = zeros(initPop,1);     % LHC sampling of initial points
            for i=1:initPop
                pop_cost(i)=CostFunction(X_all(i,:));   % evaluate initial population
            end
            clear added_FLAG;
            added_FLAG(1:initPop,1) = 'I';              % flags initial population values/points
            [maxval,maxpos] = max(pop_cost);            % find best one
            perm = randperm(length(pop_cost));          
            parents = [X_all(maxpos,:);X_all(perm(1:nPop-1),:)]; % choose parents
            k = k_est(pop_cost,X_all');                 % estimate the lipschitz constant
            best_val= []; nr_eval = [];                 
            if add_RBF > 0
                RBF_model = rbf_build(X_all, pop_cost, ker, c, pol_part, 0);   % build the RBF surrogate model
            end
            for iter = 1:MaxIt                  % simple iteration counter
                children = zeros(nPop,D); children_vals = zeros(nPop,1); LIPO_vals = zeros(nPop,1); RBF_vals = zeros(nPop,1);
                set1 = 2:nPop; 
                for i=1:nPop
                    r1 = set1(randi(nPop-1));   % choose parent 1
                    set2 = setdiff(set1,r1);
                    r2 = set2(randi(nPop-2));   % choose parent 2
                    v = parents(1,:) + F*(parents(r1,:) - parents(r2,:));   % compute v
                    v = min(v,VarMax); v = max(v,VarMin);
                    u = zeros(1,D); RN = rand(1,D);
                    for j=1:D                   % compute child position
                        if RN(j) <= CR
                           u(j) = v(j);
                        else
                           u(j) = parents(1,j);
                        end
                    end
                    if add_RBF > 0
                         RBF_vals(i)=rbf_predict(RBF_model, X_all, u);  % compute RBF surrogate value of child
                    end
                    if add_LIPO > 0 && (iter < 200 || mod(iter,5) == 0)
                        LIPO_vals(i) = fub(pop_cost,X_all',k,u');       % compute Lipschitz surrogate value of child
                    else 
                        LIPO_flag = 0;
                    end
                    children(i,:) = u;
                end
                if add_LIPO > 0 && mod(iter,f_li(iter)) == 0
                    [LIPO_v,LIPO_p] = sort(LIPO_vals,'descend');                        % find best ones (highest value) 
                    for j=1:1 %add_LIPO
                        if ~ismember(X_all, children(LIPO_p(j),:), 'rows')              % add to the whole population (if not already there)
                            X_all = [X_all;children(LIPO_p(j),:)];
                            pop_cost = [pop_cost;CostFunction(children(LIPO_p(j),:))];
                            LIPO_flag = 1;
                            added_FLAG(end+1) = 'L';    % flag point as added by L surrogate
                        else
                            LIPO_flag = 0;
                        end
                    end
                end
                if add_RBF > 0                                                          % condition on adding RBF surrogate children
                    [RBF_v,RBF_p] = sort(RBF_vals,'descend');                           % find best ones (highest value) 
                    for j=1:add_RBF                                                     
                        if ~ismember(X_all, children(RBF_p(j),:), 'rows')               % add to the whole population (if not already there)
                            X_all = [X_all;children(RBF_p(j),:)];
                            pop_cost = [pop_cost;CostFunction(children(RBF_p(j),:))];
                            RBF_flag = 1;
                            added_FLAG(end+1) = 'R';    % flag point as added by RBF surrogate
                        else
                            RBF_flag = 0;
                        end
                    end
                end
                if (add_FMINCON > 0) && mod(iter,f_lo(iter)) == 0
                    quad_vals = min(length(pop_cost),3*D);              % number of points for local surrogate model
                    [s_v,s_p] = sort(pop_cost,'descend');               % choosing the best points (highest objective value) to build local surrogate
                    X_sub = X_all(s_p(1:quad_vals),:);
                    subMax = max(X_sub); subMin = min(X_sub);           % bounds on local optimization
                    c_sub = c; ker_sub = ker; pol_part_sub = pol_part;
                    RBF_submodel = rbf_build(X_sub+0.000001*randn(size(X_sub)), s_v(1:quad_vals), ker_sub, c_sub, pol_part_sub, 0); % use polynomial part
                    func = @(x) -rbf_predict(RBF_submodel, X_sub, x);   % local objective function (-surrogate)
                    options = optimoptions("fmincon","Display","off",'maxiter',3000,'Algorithm','sqp'); %local optimization parameters (default, except for algorithm and display)
                    X0 = subMin + (subMax-subMin).*rand(1,D);
                    [xn,fval] = fmincon(func,X0,[],[],[],[],subMin,subMax,[],options); % find minimum of -surrogate (or maximum of surrogate)
                    if ~ismember(X_all, xn, 'rows')                     % add optimal point to the whole population (if not already there)
                        X_all = [X_all;xn];
                        pop_cost = [pop_cost;CostFunction(xn)];
                        added_FLAG(end+1) = 'O';        % flag point as added by local Optimization
                    end
                end
                %%
                [maxval,maxpos] = max(pop_cost);                        % find the best solution so far
                perm = randperm(length(pop_cost));
                parents = [X_all(maxpos,:);X_all(perm(1:nPop-1),:)];    % new parents
                RBF_model = rbf_build(X_all, pop_cost, ker, c, pol_part, 0);   % new RBF model
                if mod(iter,3) == 0
                    disp({-pop_cost(maxpos),length(pop_cost),added_FLAG(maxpos),trial}); % show progress (best function value, number of function evaluations, flag of best solution)
                end                                                                      % function values have a flipped sign (maximization -> minimization)
                if length(pop_cost) > MaxFe % end if over function evaluation limit
                    break;
                end
                best_val(iter) = -pop_cost(maxpos); nr_eval(iter) = length(pop_cost);
                k = k_update(pop_cost,X_all',k); % update Lipschitz constant
            end
            disp({-pop_cost(maxpos),length(pop_cost),added_FLAG(maxpos),trial});
            best_results(trial) = best_val(end);
            best_vals{trial} = best_val;
            nr_evals{trial} = nr_eval;
            com_time(trial) = toc(tStart);
        end
    end
end
end