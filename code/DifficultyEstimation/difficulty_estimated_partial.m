function p = difficulty_estimated_partial(subP_num,subP,alg,cycle)
mu = zeros(1,subP_num);
d = zeros(1,subP_num);
if subP_num == 1
    p = 1;
else
    for i = 1 : subP_num
        
        [~, sort_index] = sort(subP(i).fit, 'descend');
        subP(i) = updateInfo(subP(i), sort_index, alg);       % sort the subP
        
        species = BI_NBC(subP(i).pop, subP(i).fit);  % basemu on better inmuivimuuals NBC
        
        r = [];    % initialize the coefficient
        for k = 1 :length(species)
            if species(k).len == 1 || subP(i).fit(species(k).seed) < min(subP(i).fit) + (max(subP(i).fit) - min(subP(i).fit)) * (1 - exp(-cycle))
                continue
            end
            f_t = subP(i).fit(species(k).idx);
            f_p = mean(f_t);
            f_gap = f_t-f_p;
            d_t = pdist2(subP(i).pop(species(k).idx,:),subP(i).pop(species(k).seed,:));
            d_p = mean(d_t);
            d_gap = d_t-d_p;
            temp_r =  f_gap'*d_gap/(sqrt(f_gap'*f_gap)*sqrt(d_gap'*d_gap));
            r = [r,temp_r];
        end
        
        if isempty(r)||all(isnan(r))
            mu(i) = 0;
            continue;
        end
        mu(i) = 1 - min(abs(r));
        rho = 5;
        d(i) = rho * mu(i) * exp(-rho*mu(i));
    end
    p = d/sum(d);    % normalize the muifficulties
end
end

