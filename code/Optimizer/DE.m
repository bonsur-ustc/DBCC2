function population = DE(population,func_num, dimIndex, algRand)
[NP,D] = size(population.pop);
count = NP;   
if count < 5 
    return 
end

pop = population.pop;
fit = population.fit;

prob = get_info(func_num);
idx = ones(count,4);
F = zeros(count,2);
type1 = rand(algRand, count, 1);

% get the nearest neighbor of each individuals
dist = pdist2(pop, pop);
dist(logical(eye(NP))) = inf;
[~, nn] = min(dist, [], 2);
nn = nn(1:count); % mutate within in the nearest individuals

for i = 1 : count
    if type1(i) < 0.5 % DE/nrand/1
        idx(i, 1:2) = randperm(algRand, NP - 1, 2);
        idx(i, idx(i, :) >= i) = idx(i, idx(i, :) >= i) + 1;
        F(i, 1) = 0.5;
        F(i, 2) = 0;
    else % DE/nrand/2
        idx(i, 1:4) = randperm(algRand, NP-1, 4);
        idx(i, idx(i, :) >= i) = idx(i, idx(i, :) >= i) + 1;
        F(i, 1) = 0.5;
        F(i, 2) = 0.5;
    end
end
spring = pop(nn, :) + repmat(F(:, 1), 1, D) .* (pop(idx(:, 1), :) - pop(idx(:, 2), :)) + repmat(F(:, 2), 1, D) .* (pop(idx(:, 3), :) - pop(idx(:, 4), :));
% crossover operator
CR = 0.9;
cross = rand(algRand, count, D) < CR;
for i = 1:count
    if sum(cross(i, :)) == 0
        cross(i, randi(algRand, D)) = true;
    end
end
spring = cross .* spring + (1-cross) .* pop(1:count, :);

% check the boundary
spring = boundary_check(spring, prob.lb, prob.ub,dimIndex);
spring_val = func(spring,func_num,true,dimIndex,prob.lb);
% the selector operator

compare = fit < spring_val;
next_pop = pop;
next_pop(compare,:) = spring(compare, :);
next_fit = compare .* spring_val + (1 - compare) .* fit;

population.pop = next_pop;
population.fit = next_fit;

end