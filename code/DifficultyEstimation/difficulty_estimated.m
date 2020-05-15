function p = difficulty_estimated(subP_num,subP,pop)
%DIFFICULTY_ESTIMATED 此处显示有关此函数的摘要
%   此处显示详细说明
d = zeros(1,subP_num);

for i = 1 : subP_num
    dimIndex = subP(i).dimIndex;
    [subP(i).fit, sort_index] = sort(subP(i).fit, 'descend');
    pop(:,dimIndex) = pop(sort_index, dimIndex);
    
    TP = pop(:, dimIndex); 
    species = BI_NBC(TP, subP(i).fit);  % based on better individuals NBC
%     species = p_NBC(TP);
%     species = NBC(TP);
    
    r = [];    % initialize the coefficient
    for k = 1 :length(species)
        if species(k).len == 1
            continue;
        end
        f_t = subP(i).fit(species(k).idx);
        f_p = mean(f_t);
        f_gap = f_t-f_p;
        d_t = pdist2(pop(species(k).idx,dimIndex),pop(species(k).seed,dimIndex));
        d_p = mean(d_t);
        d_gap = d_t-d_p;
        temp_r =  f_gap'*d_gap/(sqrt(f_gap'*f_gap)*sqrt(d_gap'*d_gap)); 
        r = [r,temp_r];
    end
    d(i) = 1 - min(abs(r));
end

p = d/sum(d);    % normalize the difficulties 

end

