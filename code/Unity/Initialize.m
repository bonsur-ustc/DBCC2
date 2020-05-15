function population = Initialize(NP,prob,algRand,Alg)
%INITIALIZE 此处显示有关此函数的摘要
%   此处显示详细说明
lb = prob.lb;
ub = prob.ub;
D =  prob.dim;
population(1:NP) = struct();
if strcmp(Alg, 'PSO')            % adopte PSO algorithm
    max_v = (ub - lb)/2;
    min_v = -max_v;
    pSelf = rand(algRand, NP, D) .* (ub - lb) + lb;
    pV =  rand(algRand, NP, D) .* (max_v - min_v) + min_v;
    pBest = pSelf;

    % initialize
    for i = 1 : NP
        population(i).pop = pBest(i, :);            % pop record the pbest in PSO
        population(i).v = pV(i, :);
        population(i).cur = pSelf(i ,:);
    end
else                             % adopt other algorithms
    
    pop = rand(algRand, NP, D) .* (ub - lb) + lb;
%     pop = Initialize_population(NP, D, lb, ub, algRand);
    for i = 1 : NP
        population(i).pop = pop(i, :);            % pop record the pbest in DE & CSA
    end
    
end

