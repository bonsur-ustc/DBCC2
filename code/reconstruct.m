function solutions = reconstruct(best, subP_num, Dim, func_num, sub_prob)
%reconstruct solutions
needForSplice = 1;
for i = 1 : subP_num
    needForSplice = needForSplice * size(best(i).finalSol, 1); % The evaluations for solutons reconstructing
end


for i = 1:subP_num
    dimIndex = sub_prob{i};
    archiveX = best(i).finalSol;
    [archiveX_num,~] = size(archiveX);
    if i == 1
        solutions = zeros(archiveX_num,Dim);
        solutions(:,dimIndex) = archiveX;
    else
        [interval,~] = size(solutions); 
        solutions = repmat(solutions,archiveX_num,1);
        
        for j = 1 : archiveX_num
            row_index = (j-1)*interval+1: j*interval;
            solutions(row_index,dimIndex) = repmat(archiveX(j,:),interval,1);
        end
    end
end




