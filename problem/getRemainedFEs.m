function remain_FEs = getRemainedFEs(func_num)
global evals;
prob = get_info(func_num);
remain_FEs = prob.maxFEs - evals;
end

