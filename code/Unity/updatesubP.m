function subP = updatesubP(subP, subSpecies, index, Alg)
%UPDATESUBP 此处显示有关此函数的摘要
%   此处显示详细说明
subP.pop(index, :) = subSpecies.pop;
subP.fit(index) = subSpecies.fit;
if strcmp(Alg, 'PSO')
    subP.v(index,:) = subSpecies.v;
    subP.cur(index,:) = subSpecies.cur;
    subP.curfit(index) = subSpecies.curfit;
end
end

