function groups = Problem_seperate(func_num)
prob = get_info(func_num);
lb = prob.lb;
ub = prob.ub;

% the initial dimensional groups = {1,2,....D}
groups = 1:prob.dim;
groups = num2cell(groups);

lbFitVal= func(lb, func_num,true);

% sepearte problem
groups = fast_separating(groups, lb, lbFitVal, ub, func_num);  
end




