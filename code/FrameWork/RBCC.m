function record = RBCC(func_num, run, alg)
algRand = RandStream.create('mt19937ar','seed', run);
global evals;
global initial_flag;
initial_flag = 0; % should set the flag to 0 for each run, each function

NP = 500;                                   % population size
evals = 0;
epsilon = 0.1;
radius = 0.1;

prob = get_info(func_num);
D = prob.dim;


P = Initialize(NP,prob,algRand,alg);

variables = Problem_seperate(func_num);       % problem separating with modified RDG
subP_num = size(variables,2);

% subP records the information of subproblem 
subP = getSub(variables, subP_num, func_num, P, alg);

RemainedFEs = getRemainedFEs(func_num);
max_Gen = floor(RemainedFEs/NP);

if strcmp(alg, 'CSA')
    nc = 2;              % the number of clone is set to 2
    max_Gen = floor(RemainedFEs/(NP*nc));
end

cur_gen = 1;

%%  RBCC framework
cur = 1;
while cur_gen <= max_Gen
    disp(cur_gen);
    dimIndex = subP(cur).dimIndex;
    [~, sort_index] = sort(subP(cur).fit, 'descend');         % sort individuals by descending
    subP(cur) = updateInfo(subP(cur), sort_index, alg);                  % sort the subP
    
    species = BI_NBC(subP(cur).pop, subP(cur).fit);  % clustering with BI-NBC
    
    for i = 1 : length(species)
        % construct struct for each species
        subSpecies = creatSpe(subP(cur), species(i).idx, alg);
        % optimizing
        Lbest = creatSpe(subP(cur), species(i).seed, alg);
        if strcmp(alg, 'DE')
            subSpecies = DE(subSpecies,func_num, dimIndex, algRand);   % modified DE_nn
            % or PSO
        elseif strcmp(alg, 'PSO')
            subSpecies = PSO(subSpecies, Lbest, func_num, algRand, dimIndex);
        elseif strcmp(alg, 'CSA')
             subSpecies = CSA(subSpecies, Lbest, nc, func_num, algRand, dimIndex);
        else
            error('There is no such optimizer');
        end
        subP(cur) = updatesubP(subP(cur), subSpecies, species(i).idx, alg);
    end
    cur_gen = cur_gen + 1;
    cur = mod(cur, subP_num) + 1;
end

%%
best = struct;     % get best solution of each subproblem

for i = 1 : subP_num
    [subP(i).fit,sort_index] = sort(subP(i).fit,'descend');
    Temp_pop = subP(i).pop;
    Temp_pop = Temp_pop(sort_index,:);
    absfit = abs(subP(i).fit(1) - subP(i).fit);
    optima_set = Temp_pop((absfit<epsilon),:);
    best(i).finalSol = [];
    for k = 1 :size(optima_set,1)
        if isempty(best(i).finalSol)
            best(i).finalSol = optima_set(1,:);
            continue;
        end
        temp_dis = pdist2(optima_set(k,:),best(i).finalSol);
        if all(temp_dis > radius)
            best(i).finalSol = [best(i).finalSol;optima_set(k,:)];
        end
    end
end
 

final_solutions = reconstruct(best, subP_num, D, func_num, variables);  % reconstruct solutions


record = zeros(1, 5);

for i = 1:5
    record(i) = count_goptima(final_solutions, func_num, 10^(-i));
end
end

