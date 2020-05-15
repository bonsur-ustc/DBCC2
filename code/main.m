% This is the source code of DBCC2
% Date: 2019/5/14
% Last modified: 2019/12/19
%===============================main program==============================%
function main()
test_case = {
    @DBCC2, 'DE';
    @DBCC2, 'PSO';
    @DBCC2, 'CSA';
    
    @DBCC, 'DE';
    @DBCC, 'PSO';
    @DBCC, 'CSA';
    
    @CBCC3, 'DE';
    @CBCC3, 'PSO';
    @CBCC3, 'CSA';
    
    @RBCC, 'DE';
    @RBCC, 'PSO';
    @RBCC, 'CSA';};

test_ex = 1:size(test_case, 1);
test_func = 1:15;   % Function index

for i = 1:size(test_case,1)
    for func = test_func 
        mkdir(sprintf('./result/ALG%d',i));
        if ~ismember(func, test_func)
            continue;
        end
        delete(gcp('nocreate'));
        parpool('local',25);
        spmd(25)
            result1 = test_case{i,1}(func, labindex, test_case{i, 2}); % parallel for 25 workers
        end
        result = cat(1, result1{1:end});
        delete(gcp('nocreate'));
        parpool('local',25);
        spmd(25)
            result1 = test_case{i}(func, labindex+25, test_case{i, 2}); %parallel for 25 workers
        end
        result = [result; cat(1, result1{1:end})];  % result records the number of peaks found.
        pr = mean(result)/ get_no_goptima(func);
        stdPR = std(result./ get_no_goptima(func));
        result = [pr;stdPR; result];
        dlmwrite(sprintf('./result/ALG%d/F%d.csv',i,func), result);
    end

end
end
