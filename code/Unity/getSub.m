function subP = getSub(variables, subP_num, func_num, P, Alg)
pop = cat(1, P(:).pop);
subP = struct();
prob = get_info(func_num);
lb = prob.lb;

if strcmp(Alg, 'PSO')
    v = cat(1, P(:).v);
    cur = cat(1, P(:).cur);
    for i = 1 : subP_num
        subP(i).dimIndex = variables{i};
        subP(i).pop = pop(:, subP(i).dimIndex);
        subP(i).fit = func(subP(i).pop, func_num, true, subP(i).dimIndex, lb);
        subP(i).v = v(:, subP(i).dimIndex);
        subP(i).cur = cur(:, subP(i).dimIndex);
        subP(i).curfit = subP(i).fit;
    end
else
    for i = 1 : subP_num
        subP(i).dimIndex = variables{i};
        subP(i).pop = pop(:, subP(i).dimIndex);
        subP(i).fit = func(subP(i).pop, func_num, true, subP(i).dimIndex, lb);
    end
end
end

