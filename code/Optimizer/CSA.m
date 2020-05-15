function population = CSA(population, bestIndiv, copy, func_num, algRand, dimIndex)

[NP, D] = size(population.pop);
pop = population.pop;
fit = population.fit;
best = bestIndiv.pop;

dis = pdist2(pop, pop);
c = max(dis(:));
if NP == 1 
    c = 1;
end

alpha = 0.5 * c / sqrt(D);
prob = get_info(func_num);


pool = repmat(pop, 1, copy + 1);
pool_fit = repmat(fit, 1, copy + 1);

type = rand(algRand, NP-1,1);
% matrix operation
idx = type < 0.5 ;  % find the index of type < 0.5
idx = [0;idx];

step = 0.5*idx.*(repmat(best, NP, 1) - pool(:, 1 : D)) + alpha * (1-idx) .* randn(algRand, NP, D);
pool(:, D+1:2*D) = pool(:, D+1: 2*D) + step;
pool(:, D+1:2*D) = boundary_check(pool(:,D+1:2*D), prob.lb, prob.ub, dimIndex);

for i = 3 : copy + 1
    step = alpha .* randn(algRand, NP, D);
    pool(:, (i-1)*D+1 : i*D) = pool(:, (i-1)*D+1 : i*D) + step;
    pool(:,(i-1)*D+1:i*D) = boundary_check(pool(:,(i-1)*D+1:i*D), prob.lb, prob.ub, dimIndex);
end


for i = 2 : copy + 1
   pool_fit(:, i) = func(pool(:, (i-1)*D +1 : i*D), func_num, true, dimIndex, prob.lb); 
end

[~,max_index] = max(pool_fit, [], 2);

p_evo = zeros(NP,D);
fit_evo = zeros(NP,1);
for i = 1 : NP
    p_evo(i , :) = pool(i,(max_index(i)-1)*D +1 : max_index(i) * D);
    fit_evo(i) = pool_fit(i ,max_index(i));
end

population.pop = p_evo;
population.fit = fit_evo;

end