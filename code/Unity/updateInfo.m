function subP = updateInfo(subP, Index, Alg)
%UPDATEP 此处显示有关此函数的摘要
%   此处显示详细说明
subP.pop = subP.pop(Index, :);
subP.fit = subP.fit(Index);
if strcmp(Alg, 'PSO')
    subP.v = subP.v(Index, :);
    subP.cur = subP.cur(Index, :);
    subP.curfit = subP.curfit(Index);
end
end

