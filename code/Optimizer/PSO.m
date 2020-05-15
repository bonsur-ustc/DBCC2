function population = PSO(population, bestIndiv, func_num, algRand,dimIndex)

w = 0.729;
c1 = 2.05;
c2 = 2.05;
[NP, D] = size(population.pop);

% Set the boundary infomation 
prob = get_info(func_num);
xmin = prob.lb(dimIndex);
xmax = prob.ub(dimIndex);
vmax = (xmax - xmin)/2;
vmin = -vmax;

% Load the infomation of population 
v = population.v;
x = population.cur;
fit = population.curfit;
pbest = population.pop;             % the pbest is our mainly concern
pbestCosts = population.fit;

% seed = bestIndiv.pop;               % the seed of this population % is seed?
[~,bestIndiv] = max(pbestCosts); 
seed = pbest(bestIndiv,:);


v_tmp = w * v + w*c1 * rand(algRand, NP, D) .* (pbest-x) + w* c2 * rand(algRand, NP, D) .* (seed-x);
% Clamp veloctiy
oneForViolation = v_tmp < repmat(vmin,NP,1);
v_tmp = (1-oneForViolation).*v_tmp + oneForViolation.*repmat(vmin,NP,1);
oneForViolation = v_tmp > repmat(vmax,NP,1);
v_tmp = (1-oneForViolation).*v_tmp + oneForViolation.*repmat(vmax,NP,1);

% Update position
x_tmp = x + v_tmp;

% Reflect-Z for particles out of bounds -- S. Helwig, J. Branke, and S. Mostaghim, "Experimental Analysis of Bound Handling Techniques in Particle Swarm Optimization," IEEE TEC: 17(2), 2013, pp. 259-271

% reflect lower bound
relectionAmount = repmat(xmin,NP,1) - x_tmp;
oneForNeedReflection = relectionAmount > zeros(NP,D);
relectionAmount = rem(repmat(xmin,NP,1) - x,repmat(xmax,NP,1)-repmat(xmin,NP,1));  % amount need for reflect
relectionAmount = (1-oneForNeedReflection).*zeros(NP,D) + oneForNeedReflection.*relectionAmount;
% clampfirst
x_tmp = (1-oneForNeedReflection).*x_tmp + oneForNeedReflection.*repmat(xmin,NP,1);
% then reflect
x_tmp = x_tmp+ relectionAmount;
% set velocity for reflected particles to zero
v_tmp = (1-oneForNeedReflection).*v_tmp + oneForNeedReflection.*zeros(NP,D);

% reflect upper bound
relectionAmount = repmat(xmax,NP,1) - x_tmp;
oneForNeedReflection = relectionAmount < zeros(NP,D);
relectionAmount = rem(x - repmat(xmax,NP,1),repmat(xmax,NP,1)-repmat(xmin,NP,1)); % amount need for reflect
relectionAmount = (1-oneForNeedReflection).*zeros(NP,D) + oneForNeedReflection.*relectionAmount;
% clampfirst
x_tmp = (1-oneForNeedReflection).*x_tmp + oneForNeedReflection.*repmat(xmax,NP,1);
% then reflect
x_tmp = x_tmp - relectionAmount;
% set velocity for reflected particles to zero
v_tmp = (1-oneForNeedReflection).*v_tmp + oneForNeedReflection.*zeros(NP,D);

% evaluate, Update pbest,pbestCosts

newCosts = func(x_tmp,func_num, true, dimIndex, prob.lb);

x = x_tmp;
v = v_tmp;
improved = newCosts > pbestCosts;
pbest(improved,:) = x(improved,:);
pbestCosts(improved) = newCosts(improved);

% record the update infomation
population.v = v;
population.cur = x;
population.curfit = newCosts;
population.pop = pbest;
population.fit = pbestCosts;

end
    
